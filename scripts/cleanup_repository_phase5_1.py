#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fase 5.1 — limpeza segura do repositório.

Objetivo:
- evitar arquivos soltos na raiz;
- mover scripts geradores para scripts/generators/;
- mover corretivos pontuais para scripts/fixes/;
- arquivar READMEs temporários em docs/archive/;
- arquivar documentos Markdown antigos/sobressalentes, se solicitado;
- criar relatório de limpeza;
- preservar Fortran, LaTeX, Mathematica, papers e tese original.

Uso recomendado:

1) Dry-run, sem alterar nada:
   python scripts/cleanup_repository_phase5_1.py

2) Aplicar limpeza segura:
   python scripts/cleanup_repository_phase5_1.py --apply

3) Aplicar e arquivar backups README.before-*.md:
   python scripts/cleanup_repository_phase5_1.py --apply --archive-readme-backups

4) Aplicar e também arquivar scripts de fase/correção:
   python scripts/cleanup_repository_phase5_1.py --apply --archive-readme-backups --move-phase-scripts

Observação:
- O script não executa git add/commit automaticamente.
- Ele também não remove outputs/; apenas recomenda ignorar outputs/ via .gitignore.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


PROTECTED_PREFIXES = (
    "src/fortran_original/",
    "src/mathematica_original/",
    "thesis/original/",
    "papers/original/",
)

PROTECTED_EXTENSIONS = {
    ".f",
    ".for",
    ".f90",
    ".tex",
    ".bib",
    ".nb",
    ".eps",
    ".pdf",
}


ROOT_GENERATOR_PATTERNS = (
    "criar_*.py",
    "corrigir_*.py",
)

ROOT_SCRIPT_PATTERNS = (
    "reorganizar_*.ps1",
)

README_BACKUP_PATTERNS = (
    "README.before-*.md",
    "README.backup*.md",
    "README.old*.md",
)

TEMP_MARKDOWN_ROOT_NAMES = {
    "Texto colado.txt",
}


@dataclass
class Action:
    kind: str
    source: Path
    target: Path | None
    reason: str


def run_git(args: list[str]) -> str:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=Path.cwd(),
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            return result.stderr.strip()
        return result.stdout.strip()
    except FileNotFoundError:
        return "Git nao encontrado no PATH."


def is_repo_root(root: Path) -> bool:
    return (root / ".git").exists() and (root / "pyproject.toml").exists()


def path_posix(path: Path) -> str:
    return path.as_posix()


def is_protected(path: Path) -> bool:
    rel = path_posix(path.relative_to(Path.cwd()))
    if rel.startswith(PROTECTED_PREFIXES):
        return True
    if path.suffix.lower() in PROTECTED_EXTENSIONS and (
        "fortran_original" in path.parts
        or "mathematica_original" in path.parts
        or "thesis" in path.parts
        or "papers" in path.parts
    ):
        return True
    return False


def ensure_gitignore_entry(lines: list[str], entry: str) -> list[str]:
    stripped = [line.strip() for line in lines]
    if entry not in stripped:
        lines.append(entry)
    return lines


def update_gitignore(root: Path, apply: bool, actions: list[Action]) -> None:
    gitignore = root / ".gitignore"
    existing = gitignore.read_text(encoding="utf-8").splitlines() if gitignore.exists() else []

    desired = [
        "",
        "# Generated outputs",
        "outputs/",
        "",
        "# Python caches",
        "__pycache__/",
        "*.py[cod]",
        "*.egg-info/",
        ".pytest_cache/",
        ".mypy_cache/",
        ".ruff_cache/",
        ".venv/",
        "",
        "# OS/editor",
        ".DS_Store",
        "Thumbs.db",
        ".vscode/",
    ]

    new_lines = existing[:]
    for entry in desired:
        if entry == "":
            continue
        new_lines = ensure_gitignore_entry(new_lines, entry)

    if new_lines != existing:
        actions.append(
            Action(
                kind="update",
                source=gitignore,
                target=None,
                reason="Adicionar entradas padrao para outputs, caches Python e arquivos de editor.",
            )
        )
        if apply:
            gitignore.write_text("\n".join(new_lines).strip() + "\n", encoding="utf-8")


def collect_actions(root: Path, args: argparse.Namespace) -> list[Action]:
    actions: list[Action] = []

    if args.move_phase_scripts:
        for pattern in ROOT_GENERATOR_PATTERNS:
            for file in root.glob(pattern):
                if file.is_file():
                    if file.name.startswith("corrigir_"):
                        target = root / "scripts" / "fixes" / file.name
                        reason = "Mover script corretivo para scripts/fixes/."
                    else:
                        target = root / "scripts" / "generators" / file.name
                        reason = "Mover script gerador de fase para scripts/generators/."
                    actions.append(Action("move", file, target, reason))

        for pattern in ROOT_SCRIPT_PATTERNS:
            for file in root.glob(pattern):
                if file.is_file():
                    target = root / "scripts" / "generators" / file.name
                    actions.append(
                        Action(
                            "move",
                            file,
                            target,
                            "Mover script PowerShell de reorganizacao para scripts/generators/.",
                        )
                    )

    if args.archive_readme_backups:
        for pattern in README_BACKUP_PATTERNS:
            for file in root.glob(pattern):
                if file.is_file():
                    target = root / "docs" / "archive" / "readme-backups" / file.name
                    actions.append(
                        Action(
                            "move",
                            file,
                            target,
                            "Arquivar backup antigo de README para evitar arquivo solto na raiz.",
                        )
                    )

    if args.archive_temp_files:
        for name in TEMP_MARKDOWN_ROOT_NAMES:
            file = root / name
            if file.exists() and file.is_file():
                target = root / "docs" / "archive" / "temporary" / file.name
                actions.append(
                    Action(
                        "move",
                        file,
                        target,
                        "Arquivar arquivo temporario colado/exportado.",
                    )
                )

    update_gitignore(root, args.apply, actions)
    return actions


def apply_actions(actions: list[Action]) -> None:
    for action in actions:
        if action.kind == "move":
            assert action.target is not None
            if is_protected(action.source):
                print(f"[PROTEGIDO] ignorado: {action.source}")
                continue
            action.target.parent.mkdir(parents=True, exist_ok=True)
            if action.target.exists():
                print(f"[SKIP] destino ja existe: {action.target}")
                continue
            shutil.move(str(action.source), str(action.target))
            print(f"[MOVIDO] {action.source} -> {action.target}")

        elif action.kind == "update":
            print(f"[ATUALIZADO] {action.source}")


def build_report(root: Path, actions: list[Action], applied: bool) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    status = run_git(["status", "--short"])
    branch = run_git(["branch", "--show-current"])

    lines = [
        "# Repository Cleanup Report — Phase 5.1",
        "",
        f"- Generated at: `{now}`",
        f"- Branch: `{branch}`",
        f"- Mode: `{'apply' if applied else 'dry-run'}`",
        "",
        "## Rule",
        "",
        "This cleanup phase preserves original scientific sources:",
        "",
        "- `src/fortran_original/`",
        "- `src/mathematica_original/`",
        "- `thesis/original/`",
        "- `papers/original/`",
        "",
        "No Fortran, LaTeX, Mathematica or original thesis/paper file should be removed by this script.",
        "",
        "## Planned/applied actions",
        "",
    ]

    if actions:
        for action in actions:
            target = f" -> `{action.target}`" if action.target else ""
            lines.append(f"- `{action.kind}` `{action.source}`{target}: {action.reason}")
    else:
        lines.append("- No cleanup actions detected.")

    lines.extend(
        [
            "",
            "## Git status at report time",
            "",
            "```text",
            status or "No status output.",
            "```",
            "",
            "## Recommended validation",
            "",
            "```powershell",
            "python -m compileall src\\mms_viscoelastic",
            "python -m pytest -q",
            "python examples\\03_case12_visualization_and_ml_dataset.py",
            "git status",
            "git diff --stat",
            "```",
        ]
    )

    return "\n".join(lines) + "\n"


def write_report(root: Path, actions: list[Action], applied: bool) -> Path:
    report = root / "docs" / "maintenance" / "repository-cleanup-phase5-1.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(build_report(root, actions, applied), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Limpeza segura do repositorio da tese."
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Aplicar as acoes. Sem isso, roda em dry-run.",
    )
    parser.add_argument(
        "--move-phase-scripts",
        action="store_true",
        help="Mover scripts criar_*.py, corrigir_*.py e reorganizar_*.ps1 para scripts/.",
    )
    parser.add_argument(
        "--archive-readme-backups",
        action="store_true",
        help="Arquivar README.before-*.md e backups similares em docs/archive/.",
    )
    parser.add_argument(
        "--archive-temp-files",
        action="store_true",
        help="Arquivar arquivos temporarios conhecidos, como Texto colado.txt.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Ativa move-phase-scripts, archive-readme-backups e archive-temp-files.",
    )
    return parser.parse_args()


def main() -> None:
    root = Path.cwd()

    if not is_repo_root(root):
        raise SystemExit(
            "[ERRO] Rode este script na raiz do repositorio, onde existem .git e pyproject.toml."
        )

    args = parse_args()
    if args.all:
        args.move_phase_scripts = True
        args.archive_readme_backups = True
        args.archive_temp_files = True

    actions = collect_actions(root, args)

    print("[INFO] Fase 5.1 — limpeza segura do repositorio")
    print(f"[INFO] Modo: {'APLICAR' if args.apply else 'DRY-RUN'}")
    print(f"[INFO] Acoes detectadas: {len(actions)}")

    for action in actions:
        target = f" -> {action.target}" if action.target else ""
        print(f"  - {action.kind}: {action.source}{target}")
        print(f"    motivo: {action.reason}")

    if args.apply:
        apply_actions(actions)

    report = write_report(root, actions, args.apply)
    print(f"[OK] Relatorio escrito: {report}")

    if not args.apply:
        print("\n[INFO] Nada foi alterado. Para aplicar:")
        print("  python scripts\\cleanup_repository_phase5_1.py --apply --all")
    else:
        print("\n[OK] Limpeza aplicada. Rode testes e revise git status.")


if __name__ == "__main__":
    main()
