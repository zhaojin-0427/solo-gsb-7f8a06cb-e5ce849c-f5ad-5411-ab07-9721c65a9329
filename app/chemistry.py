"""化学式解析器。

语法（递归下降）：

    formula   := group+
    group     := element [subscript] | "(" formula ")" [subscript]
    element   := 大写字母 [小写字母]
    subscript := 1~99 的整数（不允许前导零）

省略的下标按 1 处理，圆括号可以嵌套。
"""

from __future__ import annotations


class FormulaError(ValueError):
    """化学式语法错误，position 为出错字符在式子中的下标。"""

    def __init__(self, message: str, position: int | None = None):
        super().__init__(message)
        self.position = position


def parse_formula(text: str) -> dict[str, int]:
    """把化学式解析成 {元素符号: 原子个数}。"""
    if not isinstance(text, str) or not text.strip():
        raise FormulaError("formula must be a non-empty string")
    counts, pos = _parse_group(text, 0)
    if pos != len(text):
        ch = text[pos]
        if ch == ")":
            raise FormulaError(f"unbalanced ')' at position {pos}", pos)
        raise FormulaError(f"unexpected character {ch!r} at position {pos}", pos)
    if not counts:
        raise FormulaError("formula contains no elements")
    return counts


def _parse_group(text: str, pos: int) -> tuple[dict[str, int], int]:
    """解析到 ')' 或字符串末尾为止，返回 (计数, 停止位置)。"""
    counts: dict[str, int] = {}
    while pos < len(text):
        ch = text[pos]
        if ch == ")":
            break
        if ch == "(":
            open_pos = pos
            sub, pos = _parse_group(text, pos + 1)
            if pos >= len(text) or text[pos] != ")":
                raise FormulaError(
                    f"unbalanced '(' at position {open_pos}: missing ')'", open_pos
                )
            if not sub:
                raise FormulaError(f"empty parentheses at position {open_pos}", open_pos)
            pos += 1
            multiplier, pos = _parse_subscript(text, pos)
            for symbol, n in sub.items():
                counts[symbol] = counts.get(symbol, 0) + n * multiplier
        elif "A" <= ch <= "Z":
            symbol = ch
            pos += 1
            if pos < len(text) and "a" <= text[pos] <= "z":
                symbol += text[pos]
                pos += 1
            multiplier, pos = _parse_subscript(text, pos)
            counts[symbol] = counts.get(symbol, 0) + multiplier
        else:
            raise FormulaError(
                f"unexpected character {ch!r} at position {pos}", pos
            )
    return counts, pos


def _parse_subscript(text: str, pos: int) -> tuple[int, int]:
    """解析可选下标，省略时返回 1；合法取值 1~99。"""
    start = pos
    while pos < len(text) and "0" <= text[pos] <= "9":
        pos += 1
    if pos == start:
        return 1, pos
    digits = text[start:pos]
    if len(digits) > 2 or digits[0] == "0":
        raise FormulaError(
            f"subscript {digits!r} at position {start} is out of range (1-99)", start
        )
    return int(digits), pos
