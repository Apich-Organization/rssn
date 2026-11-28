set pagination off
set debuginfod enabled off
set disable-randomization off
python
import gdb
import os

class Count(gdb.Breakpoint):
    def __init__(self, name):
        # The exact entry keeps the SysV argument registers available.
        super().__init__("*" + name, internal=True)
        self.name = name
        self.count = 0

    def stop(self):
        self.count += 1
        if capture_fse and self.name == "FSE_decompress_wksp_bmi2":
            ptr, size, capacity, log = (
                int(gdb.parse_and_eval(reg)) for reg in ("$rdx", "$rcx", "$rsi", "$r8")
            )
            header = gdb.selected_inferior().read_memory(ptr, size)
            print("FSE", capacity, log, bytes(header).hex())
        return False

capture_fse = os.getenv("RSSN_TRACE_FSE") == "1"
points = []
for symbol in (
    "FSE_decompress_wksp_bmi2",
    "HUF_decompress4X1_usingDTable_internal_fast_c_loop",
    "HUF_decompress4X2_usingDTable_internal_fast_c_loop",
    "HUF_decompress4X1_usingDTable_internal_fast_asm_loop",
    "HUF_decompress4X2_usingDTable_internal_fast_asm_loop",
    "ZSTD_decompressSequences",
    "ZSTD_decompressSequencesLong",
    "ZSTD_decompressSequencesSplitLitBuffer",
):
    try:
        points.append(Count(symbol))
    except gdb.error as error:
        print(symbol, "unavailable:", error)
end
run
python
for point in points:
    print(point.name, point.count)
end
