import logging
from rich.console import Console

# Create a shared console instance
console = Console()

def log_step(tag, message, color="green"):
    """
    Prints a beautifully formatted, colored message to the terminal using rich,
    and a plain-text message to the standard python logger.
    
    Args:
        tag (str): The tag to display (e.g., "WORKFLOW"). If None, the tag is omitted.
        message (str): The main log message.
        color (str): The rich color string (e.g., "green", "cyan", "yellow").
    """
    # 1. Print the pretty version to the terminal
    if tag:
        console.print(f"[{color}][bold]\\[{tag}][/bold] {message}[/{color}]")
        logging.info(f"[{tag}] {message}")
    else:
        if color == "white":
            console.print(message)
        else:
            console.print(f"[{color}]{message}[/{color}]")
        logging.info(message)
