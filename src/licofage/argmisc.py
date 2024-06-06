from argparse import HelpFormatter
import textwrap


def wrap_paragraphs(text: str, width: int, indent: str):
    """
    https://stackoverflow.com/a/75344754
    Wrapper around `textwrap.wrap()` which keeps newlines in the input string
    intact.
    """
    lines = list[str]()

    for i in text.splitlines():
        paragraph_lines = textwrap.wrap(
            i, width, initial_indent=indent, subsequent_indent=indent
        )

        # `textwrap.wrap()` will return an empty list when passed an empty
        # string (which happens when there are two consecutive line breaks in
        # the input string). This would lead to those line breaks being
        # collapsed into a single line break, effectively removing empty lines
        # from the input. Thus, we add an empty line in that case.
        lines.extend(paragraph_lines or [""])

    return lines


class Formatter(HelpFormatter):
    def _split_lines(self, text, width):
        return wrap_paragraphs(text, width, "")

    def _fill_text(self, text, width, indent):
        return "\n".join(wrap_paragraphs(text, width, indent))
