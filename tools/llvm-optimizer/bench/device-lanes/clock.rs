// The public clock_gettime symbol uses the legacy timespec on 32-bit ARM
// and x86. The time64 interfaces on these ABIs use different symbols.
#[cfg(all(target_os = "linux", any(target_arch = "arm", target_arch = "x86")))]
type ClockField = i32;
#[cfg(any(
    all(
        target_os = "linux",
        not(any(target_arch = "arm", target_arch = "x86"))
    ),
    all(target_os = "freebsd", target_pointer_width = "64"),
))]
type ClockField = i64;

#[cfg(any(
    target_os = "linux",
    all(target_os = "freebsd", target_pointer_width = "64")
))]
#[repr(C)]
struct Timespec {
    sec: ClockField,
    ns: ClockField,
}
#[cfg(any(
    target_os = "linux",
    all(target_os = "freebsd", target_pointer_width = "64")
))]
unsafe extern "C" {
    fn clock_gettime(
        clock: i32,
        ts: *mut Timespec,
    ) -> i32;
}
#[cfg(any(
    target_os = "linux",
    all(target_os = "freebsd", target_pointer_width = "64")
))]
pub(super) fn now() -> u64 {
    let mut ts = Timespec { sec: 0, ns: 0 };
    let thread_cpu_clock = if cfg!(target_os = "freebsd") {
        14
    } else {
        3
    };
    assert_eq!(unsafe { clock_gettime(thread_cpu_clock, &mut ts) }, 0);
    ts.sec as u64 * 1_000_000_000 + ts.ns as u64
}
#[cfg(windows)]
pub(super) fn now() -> u64 {
    #[repr(C)]
    #[derive(Default)]
    struct FileTime {
        low: u32,
        high: u32,
    }

    #[link(name = "kernel32", kind = "raw-dylib")]
    unsafe extern "system" {
        fn GetCurrentThread() -> isize;
        fn GetThreadTimes(
            thread: isize,
            created: *mut FileTime,
            exited: *mut FileTime,
            kernel: *mut FileTime,
            user: *mut FileTime,
        ) -> i32;
    }

    let mut created = FileTime::default();
    let mut exited = FileTime::default();
    let mut kernel = FileTime::default();
    let mut user = FileTime::default();
    assert_ne!(
        unsafe {
            GetThreadTimes(
                GetCurrentThread(),
                &mut created,
                &mut exited,
                &mut kernel,
                &mut user,
            )
        },
        0
    );
    let ticks = |value: FileTime| (u64::from(value.high) << 32) | u64::from(value.low);
    100 * (ticks(kernel) + ticks(user))
}
