"""Archive-independent validation of public numeric arguments."""
from numbers import Integral, Real
import math


def object_id(value):
    text = str(value)
    digits = text[1:] if text.startswith('-') else text
    if not digits or not digits.isascii() or not digits.isdecimal():
        raise ValueError('object_id must be a decimal integer with an optional leading minus sign')
    number = int(text)
    if not -(2**63) <= number < 2**63:
        raise ValueError('object_id is outside the signed 64-bit range')
    return str(number)


def positive_scalar(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a positive finite scalar')
    return value


def row_limit(value):
    if value is not None and (isinstance(value, bool) or not isinstance(value, Integral) or value < -1):
        raise ValueError('limit must be -1, a nonnegative integer, or None')
    return value
