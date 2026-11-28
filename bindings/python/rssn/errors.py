"""Statuses and the exception the package raises.

`rssn.h` reports failures as an `RssnStatus` code plus a thread-local message.
The wrapper turns that pair into an :class:`RssnError`, so the code is available
as :attr:`RssnError.status` and the message as the exception text.
"""

from __future__ import annotations

from enum import IntEnum

__all__ = ["RssnError", "Status"]


class Status(IntEnum):
    """`RssnStatus` from `rssn.h`."""

    OK = 0
    NULL_ARGUMENT = 1
    INVALID_UTF8 = 2
    PARSE = 3
    UNKNOWN_OPERATOR = 4
    NOT_NUMERIC = 5
    INCONSISTENT = 6
    RULE = 7
    NO_SUCH_TERM = 8
    NO_VALUE = 9
    SIMULATION = 10
    PANIC = 11
    DIMENSION = 12
    INVALID_ARGUMENT = 13
    CANCELLED = 14


class RssnError(Exception):
    """A failed rssn call.

    Attributes:
        status: the :class:`Status` code of the failure.
        message: the message rssn reported, verbatim.
    """

    def __init__(self, status: Status, message: str) -> None:
        self.status = status
        self.message = message
        super().__init__(f"{status.name.lower()}: {message}" if message else status.name.lower())

    def __reduce__(self):  # keep pickling working
        return (type(self), (self.status, self.message))


_STATUS_BY_CODE = {int(status): status for status in Status}


def error_for_status(code: int, message: str | None = None) -> RssnError:
    """The exception for a numeric status code, reading the message when omitted."""
    from . import _ffi

    status = _STATUS_BY_CODE.get(int(code), Status.INVALID_ARGUMENT)
    if message is None:
        message = _ffi.last_error()
    return RssnError(status, message)


def check(code: int) -> None:
    """Raise :class:`RssnError` unless `code` is `OK`."""
    if int(code) != 0:
        raise error_for_status(code)
