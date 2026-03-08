import json
import ctypes
import os

lib_path = "./target/release/librssn.so"
if not os.path.exists(lib_path):
    print("Building librssn.so...")
    os.system("cargo build --all-features --release")

lib = ctypes.CDLL(lib_path)

# Set argtypes and restypes for all FFI functions
lib.rssn_jit_create.restype = ctypes.c_void_p
lib.rssn_jit_create.argtypes = []

lib.rssn_jit_configure_sandbox_json.restype = ctypes.c_char_p
lib.rssn_jit_configure_sandbox_json.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

lib.rssn_jit_compile_json.restype = ctypes.c_char_p
lib.rssn_jit_compile_json.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

lib.rssn_jit_execute.restype = ctypes.c_double
lib.rssn_jit_execute.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

engine = lib.rssn_jit_create()

# Configure the sandbox
conf_payload = {
    "memory_regions": [],
    "allowed_calls": []
}
json_conf = json.dumps(conf_payload).encode('utf-8')
lib.rssn_jit_configure_sandbox_json(engine, json_conf)

malicious_payload = {
    "instructions": [
        {"ImmI": 0},
        {"Load": "I64"},
        {"Return": None}
    ]
}

json_data = json.dumps(malicious_payload).encode('utf-8')
func_ptr_str = lib.rssn_jit_compile_json(engine, json_data)
if not func_ptr_str:
    print("Compilation failed")
    exit(1)

result_dict = json.loads(func_ptr_str.decode('utf-8'))

if "ok" in result_dict and result_dict["ok"] is not None:
    func_ptr = result_dict["ok"]
    print(func_ptr)
else:
    print("Compilation error:", result_dict.get("err"))
    exit(1)

print("Attempting to execute malicious JIT code...")
try:
    result = lib.rssn_jit_execute(engine, func_ptr)
    print(f"JIT execute results: {result}")
except Exception as e:
    print("Failed to execute:", e)

print("done")
