from argparse import (
    ArgumentParser,
    BooleanOptionalAction,
    Action,
    ArgumentError,
)
from pathlib import Path as P
from .api import Study
from .argmisc import Formatter
from . import __version__


class valid_dir(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        path = P(values[0])
        if not path.is_dir():
            raise ArgumentError(self, f"{values[0]} is not an existing directory")
        setattr(namespace, self.dest, path)


class valid_vector(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            U = list(map(int, values[0].strip().split(",")))
            setattr(namespace, self.dest, U)
        except:
            raise ArgumentError(self, f"{values[0]} is not a valid vector")


class valid_equation(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            fname = P(values[1])
        except:
            raise ArgumentError(self, f"{values[1]} is not a valid filename")
        try:
            vec, cst = values[0].strip().split("=")
            b = int(cst.strip())
            aa = list(map(int, vec.strip().split(",")))
            setattr(namespace, self.dest, (aa, b, fname))
        except:
            raise ArgumentError(self, f"{values[0]} is not a valid vector")


def main():
    parser = ArgumentParser(
        prog="licofage",
        formatter_class=Formatter,
        description="Generate deterministic finite automata associated to linear equations expressed in the numeration system associated to a given substitution.",
        epilog="""You might want to try one of the following substitutions:
     - Fibonacci "a->ab, b->a" (X^2-X-1)
     - Tribonacci "a->ab, b->ca, c->a" (X^3-X^2-X-1)
     - Some X*Pisot "a->abba, b->ab" (X^2-3X)
     - Some larger example "a->ac, b->abc, c->bcc" (X^3-4X^2+4X-2)
    
    For Walnut usage, a typical invocation is:
    $ licofage -W xpisot -D ~/Walnutdir/ "a->abba, b->ab"
    """,
    )

    parser.add_argument(
        "-V", "--version", action="version", version="%(prog)s " + __version__
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="display more messages"
    )
    parser.add_argument(
        "-S",
        "--stats",
        action="store_true",
        help="display statistics about computations",
    )
    parser.add_argument(
        "-W",
        "--Walnut",
        nargs=1,
        metavar="NAME",
        help="generate and organize output for Walnut usage under numeration system name NAME.",
    )
    parser.add_argument(
        "-D",
        "--outdir",
        action=valid_dir,
        default=P(""),
        nargs=1,
        metavar="DIRECTORY",
        help="generate output relatively to existing directory DIRECTORY",
    )
    parser.add_argument(
        "--simplify",
        action=BooleanOptionalAction,
        default=True,
        help="simplify the substitution before study (default: True)",
    )
    parser.add_argument(
        "--minimize",
        action=BooleanOptionalAction,
        default=True,
        help="minimize the automata before output (default: True)",
    )
    parser.add_argument(
        "--revector",
        action=BooleanOptionalAction,
        default=False,
        help="try to (re)associate vectors to the minimal automaton (default: False)",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-M",
        "--MSD",
        action="store_true",
        default=True,
        help="output automata for MSD numeration (default: MSD)",
    )
    group.add_argument(
        "-L",
        "--LSD",
        action="store_false",
        dest="MSD",
        help="output automata for LSD numeration",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["dot", "Walnut"],
        default="dot",
        help="select output automata format",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-U",
        "--vector",
        nargs=1,
        action=valid_vector,
        metavar='"x1, x2, ..., xk"',
        help="use this vector to generate the bounds",
    )
    group.add_argument(
        "-B",
        "--bound",
        nargs=1,
        action=valid_vector,
        metavar='"x1, x2, ..., xk"',
        help="use this vector as manual bounds",
    )
    parser.add_argument("subst", help="substitution considered")
    parser.add_argument(
        "-d",
        "--dfao",
        nargs="?",
        metavar="filename",
        const=P("morphism.txt"),
        type=P,
        help="generate the substitution automaton with output",
    )
    parser.add_argument(
        "-P",
        "--parent",
        nargs="?",
        metavar="filename",
        const=P("parent.txt"),
        type=P,
        help="generate the parent substitution automaton with output",
    )

    parser.add_argument(
        "-n",
        "--numsys",
        nargs="?",
        metavar="filename",
        const=P("numeration.txt"),
        type=P,
        help="generate the numeration system automaton",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-a",
        "--addition",
        nargs="?",
        metavar="filename",
        const=P("addition.txt"),
        type=P,
        help='generate the additon automaton (shortcut for --linear "1, 1, -1 = 0")',
    )
    group.add_argument(
        "-l",
        "--linear",
        default=None,
        action=valid_equation,
        nargs=2,
        metavar=('"a1, a2, ..., ak = b"', "filename"),
        help="generate the automaton of the given linear equation \sum_i a_i x_i = b. (this option can be used multiple times)",
    )
    parser.add_argument(
        "-s",
        "--seq",
        nargs="?",
        metavar="filename",
        const=P("sequence.dot"),
        type=P,
        help="generate the sequence automaton of the linear equation/addition",
    )
    parser.add_argument(
        "-p",
        "--poly",
        nargs="?",
        metavar="filename",
        const=P("poly.dot"),
        type=P,
        help="generate the polynomial automaton of the linear equation/addition",
    )

    parser.add_argument(
        "-c",
        "--check",
        nargs="?",
        metavar="filename",
        const=P("check.txt"),
        type=P,
        help="generate Walnut addition verification script",
    )

    def genparentdir(fname):
        parent = fname.parent
        if not parent.is_dir():
            parent.mkdir(parents=True)
        return fname

    name = None
    args = parser.parse_args()
    if args.Walnut is not None:
        args.format = "Walnut"
        name = args.Walnut[0].lower()
        endian = "msd_" if args.MSD else "lsd_"
        args.dfao = P("Word Automata Library") / P(f"{name.title()}.txt")
        args.parent = P("Word Automata Library") / P(f"{name.title()}Parent.txt")
        args.numsys = P("Custom Bases") / P(f"{endian}{name}.txt")
        if args.linear is None:
            args.addition = P("Custom Bases") / P(f"{endian}{name}_addition.txt")
        args.check = P("Command Files") / P(f"check_{name}.txt")

    try:
        if args.check is not None and name is None:
            raise ValueError("--check only available with --Walnut!")
        stud = Study(args.subst, args.verbose, args.stats, args.simplify)
        if args.dfao is not None:
            stud.gendfao(genparentdir(args.outdir / args.dfao), args.format)
        if args.parent is not None:
            stud.genparentdfao(genparentdir(args.outdir / args.parent), args.format)
        if args.numsys is not None:
            stud.gennumsys(
                genparentdir(args.outdir / args.numsys),
                args.MSD,
                args.format,
                args.minimize,
            )
        if args.addition is not None:
            args.linear = ([1, 1, -1], 0, args.addition)
        if args.linear is not None:
            (aa, b, fname) = args.linear
            stud.genlinear(
                genparentdir(args.outdir / fname),
                (aa, b),
                args.MSD,
                args.format,
                args.minimize,
                genparentdir(args.outdir / args.seq) if args.seq is not None else None,
                genparentdir(args.outdir / args.poly)
                if args.poly is not None
                else None,
                args.vector,
                args.bound,
                args.revector,
            )
        if args.check is not None:
            stud.gencheck(genparentdir(args.outdir / args.check), args.MSD, name)
        if args.verbose:
            print("Done.")
    except ValueError as e:
        print(f"*** {str(e)}. Sorry!")
    except AssertionError as e:
        print(f"### {str(e)}. Sorry!")
