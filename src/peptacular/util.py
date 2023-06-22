from typing import List, Any


def check_parentheses(text):
    """
    Check if parentheses in the given text are balanced or not.

    The function uses a stack data structure to ensure each opening parenthesis '('
    has a corresponding closing parenthesis ')' and vice versa.

    Args:
        text (str): The input string to check for balanced parentheses.

    Returns:
        bool: True if parentheses are balanced, False otherwise.
    """
    stack = []
    for i, char in enumerate(text):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if not stack:
                return False
            else:
                stack.pop()
    return len(stack) == 0


def flatten(nested_list: List[List[Any]]) -> List[Any]:
    """
    returns a flattened versions of the nested list
    """
    return [e for sublist in nested_list for e in sublist]
