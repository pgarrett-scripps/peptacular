def check_parentheses(text):
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