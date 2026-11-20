"""本地门控. 在真实 MRC 文件头上检查 CA2 输入单位, 并验证后继工具失败时仍保留已生成受体."""

from pathlib import Path
from subprocess import CalledProcessError, CompletedProcess
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import mrcfile
import numpy as np

from ligandseek_prepare.pipeline import prepare_sample
from ligandseek_prepare.receptor import prepare_receptor


class PrepareTests(unittest.TestCase):
    """核对 CryoAtom2 实际接收的网格与样本文件的发布时机."""

    def test_standard_mrc_adapts_only_copy_given_to_cryoatom(self):
        """用非零原点和非均匀步长验证坐标转换, 原始 MRC 文件在调用后保持原值."""
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            image, fasta, out = root / "map.mrc", root / "protein.fa", root / "prepared"
            out.mkdir()
            fasta.write_text(">seq\nAGS\n", encoding="utf-8")
            values = np.arange(48, dtype=np.float32).reshape(2, 4, 6)
            with mrcfile.new(image) as handle:
                handle.set_data(values)
                handle.voxel_size = (2.0, 3.0, 4.0)
                handle.header.origin = (10.0, 15.0, 24.0)
                handle.header.nxstart, handle.header.nystart, handle.header.nzstart = (7, 6, 5)

            def official_call(command, **kwargs):
                """读取传入 CryoAtom2 的 MRC, 然后按其真实目录约定产生一份测试 CIF."""
                # 输入给 CA2 的 ORIGIN 单位是网格步长, 叠加的 NSTART 应归零.
                native = Path(command[command.index("--map-path") + 1])
                with mrcfile.open(native) as handle:
                    self.assertTrue(np.array_equal(handle.data, values))
                    self.assertEqual(tuple(handle.header.origin.tolist()), (5.0, 5.0, 6.0))
                    self.assertEqual((handle.header.nxstart, handle.header.nystart, handle.header.nzstart), (0, 0, 0))
                self.assertEqual(command[command.index("--protein-sequence-path") + 1], str(fasta))
                directory = Path(command[command.index("--output-dir") + 1])
                directory.mkdir()
                (directory / "cryoatom2.cif").write_text("data_ca2\n")
                return CompletedProcess(command, 0)

            with patch("ligandseek_prepare.receptor.subprocess.run", side_effect=official_call), patch(
                "ligandseek_prepare.receptor.subprocess.check_output", return_value="2.1.1\n"
            ):
                record = prepare_receptor({"protein_fasta": fasta.name}, image, root, out)
            self.assertEqual(record["mode"], "with_sequences")
            self.assertEqual((out / "receptor.cif").read_text(), "data_ca2\n")
            with mrcfile.open(image) as handle:
                np.testing.assert_array_equal(handle.header.origin.tolist(), [10.0, 15.0, 24.0])
                self.assertEqual(int(handle.header.nxstart), 7)

    def test_failed_worker_keeps_receptor_and_does_not_publish_sample(self):
        """准备进程中断时, 已写出的受体可供核查, 最终样本目录只在全部成功后出现."""
        with TemporaryDirectory() as tmp:
            root = Path(tmp)

            def stop_after_receptor(command, **kwargs):
                """模拟外部进程在写出受体之后返回非零退出码."""
                destination = Path(command[-1])
                (destination / "receptor.cif").write_text("data_retained\n")
                CompletedProcess(command, 1).check_returncode()

            with patch("ligandseek_prepare.pipeline.subprocess.run", side_effect=stop_after_receptor):
                with self.assertRaises(CalledProcessError):
                    prepare_sample(root / "raw.json", root / "models", root / "prepared")
            self.assertFalse((root / "prepared").exists())
            work, = root.glob("prepared.unfinished-*")
            self.assertEqual((work / "receptor.cif").read_text(), "data_retained\n")
            self.assertFalse((work / "sample.json").exists())
