from hashlib import sha256


def upload_signature(*uploads):
    """Identify uploaded contents without changing the files' read positions."""
    return tuple(sha256(upload.getvalue()).hexdigest() for upload in uploads)
