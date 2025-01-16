import os, time, sys

def compute_hash(filepath):
    result = 0
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(128 * 1024), b""):
            result = hash(chunk + abs(result).to_bytes(8, 'big'))
    return result


hashes = {}
timestamps = {}

while True:
    rebuild = False
    for root, dirs, files in os.walk(os.getcwd()):
        for file in files:
            ext = file.split('.')
            ext = ext[-1] if len(ext) > 0 else ''
            if ext in ['cc', 'h', 'json', 'frag', 'vert', 'BUILD', 'bazel', 'comp'] or file in ['BUILD', 'WORKSPACE']:
                filepath = os.path.join(root, file)
                timestamp = os.path.getmtime(filepath)
                if filepath not in timestamps or timestamps[filepath] != timestamp:
                    hash_code = compute_hash(filepath)
                    if filepath not in hashes or hashes[filepath] != hash_code:
                        hashes[filepath] = hash_code
                        rebuild = True

    if rebuild:
        os.system('clear')
        os.system(" ".join(sys.argv[1:]))

    time.sleep(1)
