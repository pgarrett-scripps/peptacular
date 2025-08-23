from enum import StrEnum


class ProformaChar(StrEnum):
    MOD_START = "["
    MOD_END = "]"
    STATIC_MOD_START = "<"
    STATIC_MOD_END = ">"
    INTERVAL_START = "("
    INTERVAL_END = ")"
    UNKNOWN = "?"
    STATIC_SEP = "@"
    LABILE_MOD_START = "{"
    LABILE_MOD_END = "}"
    TERM_MOD = "-"
    CHIMERIC = '+'
    CONNECTED = '/'
    PLUS = '+'
    MINUS = '-'
    MULTI = '^'
    CHARGE_SEP = '/'