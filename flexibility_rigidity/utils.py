def print_progress(current, total, message="Processing"):
    """
    Print a progress bar for tracking progress.
    """
    progress = (current / total) * 100
    print(f"{message}: {progress:.2f}%", end="\r")
