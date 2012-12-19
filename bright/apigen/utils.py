"""Helper functions for bright API generation."""

def indent(s, n):
    """Indents all lines in the string or list s by n spaces."""
    spaces = " " * n
    lines = s.splitlines() if isinstance(s, basestring) else s
    return '\n'.join([spaces + l for l in lines])
