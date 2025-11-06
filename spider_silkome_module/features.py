import os
from pathlib import Path
import subprocess

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import PROCESSED_DATA_DIR

app = typer.Typer()


def run_shell_command_with_check(command: str, output_file: str, force: bool = False) -> bool:
    """
    Execute shell commands to automatically check if the output file exists.

    If the output file already exists, skip execution; if it does not exist, execute the command.

    Parameters
    ----------
    command : str
        Shell command to execute
    output_file : str
        Output file path of the command
    force : bool, optional
        Whether to force the execution of the command (even if the output file already exists), defaults to False.

    Returns
    -------
    bool
        True indicates that the command is executed, False indicates that the command is skipped.

    Examples
    --------
    >>> cmd = "grep -E 'gene' input.gff > output.gff"
    >>> run_shell_command_with_check(cmd, "output.gff")

    Raises
    ------
    subprocess.CalledProcessError
        If the command execution fails
    """
    # 检查输出文件是否存在
    if os.path.exists(output_file) and not force:
        logger.info(f"The output file already exists, skipping execution: {output_file}")
        return False

    # 执行命令
    logger.info(f"Execute command: {command}")
    try:
        subprocess.run(command, shell=True)
        logger.success(f"Command executed successfully, output file: {output_file}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command execution failed: {e}")
        raise


@app.command()
def main(
    # ---- REPLACE DEFAULT PATHS AS APPROPRIATE ----
    input_path: Path = PROCESSED_DATA_DIR / "dataset.csv",
    output_path: Path = PROCESSED_DATA_DIR / "features.csv",
    # -----------------------------------------
):
    # ---- REPLACE THIS WITH YOUR OWN CODE ----
    logger.info("Generating features from dataset...")
    for i in tqdm(range(10), total=10):
        if i == 5:
            logger.info("Something happened for iteration 5.")
    logger.success("Features generation complete.")
    # -----------------------------------------


if __name__ == "__main__":
    app()
