"""Publish HTML reports to Cloudflare Pages with an auto-generated index page.

First-time setup (one shot):
    npx wrangler@latest login

Default behavior: publishes every top-level ``*.html`` under ``reports/`` and
renders a landing ``index.html`` linking to each report. Idempotent: the
Pages project is auto-created on first run and reused on subsequent runs, so
the public URL stays stable for sharing.

Usage:
    pixi run python scripts/publish_report.py
    pixi run python scripts/publish_report.py reports/agent_vs_manual_report.html
    pixi run python scripts/publish_report.py reports/foo.html reports/bar.html
    pixi run python scripts/publish_report.py --project-name=my-name reports/*.html

Output URL: https://<project-name>.pages.dev/
"""

from __future__ import annotations

from datetime import datetime
import html
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Optional

from loguru import logger
import typer

app = typer.Typer(add_completion=False)

# Cloudflare project naming rules: lowercase alnum + hyphens, 1-58 chars,
# cannot start or end with a hyphen.
_PROJECT_NAME_RE = re.compile(r"^[a-z0-9](?:[a-z0-9-]{0,56}[a-z0-9])?$")
_TITLE_RE = re.compile(r"<title>([^<]+)</title>", re.IGNORECASE)

_WRANGLER_CMD = ["npx", "--yes", "wrangler@latest"]
_DEFAULT_PROJECT = "spider-silkome-reports"
_REPORTS_DIR = Path("reports")


# ─── Index page ──────────────────────────────────────────────────────────

_INDEX_TEMPLATE = """<!doctype html>
<html lang="zh-CN"><head>
<meta charset="utf-8">
<title>Spider Silkome Reports</title>
<style>
  * {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, "PingFang SC", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
         margin: 0; padding: 0; color: #222; background: #fafafa; line-height: 1.5; }}
  header {{ background: linear-gradient(135deg, #2c3e50 0%, #4a6572 100%); color: white;
           padding: 32px 40px; }}
  header h1 {{ margin: 0 0 6px 0; font-size: 26px; font-weight: 600; }}
  header .subtitle {{ opacity: 0.85; font-size: 13px; }}
  main {{ max-width: 900px; margin: 0 auto; padding: 32px 40px 60px; }}
  .card {{ display: block; background: white; border-radius: 8px; padding: 18px 22px;
          margin: 14px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.08); text-decoration: none;
          color: inherit; transition: transform 0.1s, box-shadow 0.1s; border-left: 4px solid #4477AA; }}
  .card:hover {{ transform: translateX(4px); box-shadow: 0 3px 8px rgba(0,0,0,0.12); }}
  .card .title {{ font-size: 16px; font-weight: 600; color: #2c3e50; margin: 0; }}
  .card .meta {{ font-size: 12px; color: #888; margin-top: 4px; }}
  .card .arrow {{ float: right; color: #4477AA; font-weight: 600; }}
  .empty {{ text-align: center; color: #888; padding: 40px; font-style: italic; }}
  footer {{ text-align: center; color: #aaa; font-size: 12px; padding: 20px; }}
</style>
</head><body>
<header>
  <h1>Spider Silkome Reports</h1>
  <div class="subtitle">{n_reports} report{plural} · last updated {now}</div>
</header>
<main>
{cards}
</main>
<footer>spider_silkome | published via Cloudflare Pages</footer>
</body></html>
"""

_CARD_TEMPLATE = (
    '<a class="card" href="{href}">'
    '<span class="arrow">Open →</span>'
    '<div class="title">{title}</div>'
    '<div class="meta">Updated: {mtime}　·　{size}</div>'
    "</a>"
)


def _extract_title(path: Path) -> str:
    """Best-effort parse of <title>...</title>; fall back to filename stem."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return path.stem
    m = _TITLE_RE.search(text)
    return m.group(1).strip() if m else path.stem


def _format_size(num_bytes: int) -> str:
    if num_bytes >= 1024 * 1024:
        return f"{num_bytes / 1024 / 1024:.1f} MB"
    return f"{num_bytes / 1024:.1f} KB"


def _render_card(path: Path) -> str:
    stat = path.stat()
    return _CARD_TEMPLATE.format(
        href=html.escape(path.name),
        title=html.escape(_extract_title(path)),
        mtime=datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M"),
        size=_format_size(stat.st_size),
    )


def _generate_index(deploy_dir: Path) -> int:
    """Write index.html into ``deploy_dir`` listing all *.html siblings."""
    htmls = sorted(p for p in deploy_dir.glob("*.html") if p.name != "index.html")
    cards_html = (
        "\n".join(_render_card(p) for p in htmls)
        if htmls
        else '<div class="empty">No reports published yet.</div>'
    )
    page = _INDEX_TEMPLATE.format(
        n_reports=len(htmls),
        plural="s" if len(htmls) != 1 else "",
        now=datetime.now().strftime("%Y-%m-%d %H:%M"),
        cards=cards_html,
    )
    (deploy_dir / "index.html").write_text(page, encoding="utf-8")
    return len(htmls)


# ─── Wrangler operations ─────────────────────────────────────────────────


def _run_wrangler(*args: str, capture: bool = False) -> str:
    """Run ``wrangler`` with given args. Returns stdout+stderr if capture else ''."""
    cmd = [*_WRANGLER_CMD, *args]
    if capture:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return proc.stdout + proc.stderr
    subprocess.run(cmd, check=True)
    return ""


def _check_auth() -> None:
    if "You are logged in" not in _run_wrangler("whoami", capture=True):
        logger.error("Not logged in to Cloudflare. Run once: npx wrangler@latest login")
        raise typer.Exit(1)


def _ensure_project(name: str) -> None:
    """Reuse the Pages project if it exists; otherwise create it.

    ``pages deploy --project-name`` does NOT auto-create the project — it
    errors with "Project not found". So we probe ``pages project list`` and
    create on demand, matching the name as a standalone token to avoid
    substring matches.
    """
    output = _run_wrangler("pages", "project", "list", capture=True)
    if re.search(rf"(?<!\S){re.escape(name)}(?!\S)", output):
        logger.info(f"Project '{name}' found, reusing.")
        return
    logger.info(f"Project '{name}' not found, creating...")
    _run_wrangler("pages", "project", "create", name, "--production-branch=main")


# ─── CLI entry ───────────────────────────────────────────────────────────


@app.command(help=__doc__)
def main(
    files: Optional[list[Path]] = typer.Argument(
        None,
        help="HTML files to publish (default: every reports/*.html).",
        show_default=False,
    ),
    project_name: str = typer.Option(
        _DEFAULT_PROJECT,
        "--project-name",
        help="Cloudflare Pages project name. Lowercase alnum + hyphens, 1-58 chars.",
    ),
):
    # 1. Validate project name format upfront.
    if not _PROJECT_NAME_RE.match(project_name):
        logger.error(
            f"Invalid project name '{project_name}'. Must be lowercase alphanumeric "
            "with hyphens (1-58 chars, cannot start or end with a hyphen)."
        )
        raise typer.Exit(1)

    # 2. Resolve & validate input files.
    if not files:
        files = sorted(_REPORTS_DIR.glob("*.html"))
        if not files:
            logger.error(f"No HTML files found in {_REPORTS_DIR}/.")
            logger.info("Generate one first, e.g.:")
            logger.info("  pixi run python -m spider_silkome_module.build_comparison_report")
            raise typer.Exit(1)
    for f in files:
        if not f.is_file():
            logger.error(f"File not found: {f}")
            raise typer.Exit(1)

    # 3. Auth & project pre-checks (fail fast before any staging).
    logger.info("[1/4] Checking Cloudflare auth...")
    _check_auth()

    logger.info(f"[2/4] Ensuring Pages project '{project_name}' exists...")
    _ensure_project(project_name)

    # 4. Stage to a clean temp dir, generate index, deploy.
    with tempfile.TemporaryDirectory(prefix="cf-publish-") as tmpdir:
        dist = Path(tmpdir)
        logger.info(f"[3/4] Staging {len(files)} HTML file(s)...")
        for f in files:
            shutil.copy(f, dist / f.name)
            logger.info(f"      + {f.name}")
        n = _generate_index(dist)
        logger.info(f"      + index.html ({n} card{'s' if n != 1 else ''})")

        logger.info(f"[4/4] Deploying to project '{project_name}'...")
        try:
            _run_wrangler(
                "pages", "deploy", str(dist),
                f"--project-name={project_name}",
                "--commit-dirty=true",
            )
        except subprocess.CalledProcessError as exc:
            logger.error(f"Deployment failed (exit code {exc.returncode}).")
            raise typer.Exit(1) from None

    logger.success(f"Deployed: https://{project_name}.pages.dev/")


if __name__ == "__main__":
    app()
