import sys
import unittest
from io import BytesIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from upload_signature import upload_signature


class UploadSignatureTests(unittest.TestCase):
    def test_same_size_replacement_invalidates_signature(self):
        self.assertNotEqual(upload_signature(BytesIO(b"1,2")),
                            upload_signature(BytesIO(b"3,4")))

    def test_read_position_does_not_change_signature(self):
        upload = BytesIO(b"sample,data")
        expected = upload_signature(upload)
        upload.read(3)
        self.assertEqual(expected, upload_signature(upload))
        self.assertEqual(upload.tell(), 3)


if __name__ == "__main__":
    unittest.main()
