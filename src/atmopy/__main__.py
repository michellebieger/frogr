"""Command-line interface."""
import click


@click.command()
@click.version_option()
def main() -> None:
    """ATMOpy."""


if __name__ == "__main__":
    main(prog_name="atmopy")  # pragma: no cover
