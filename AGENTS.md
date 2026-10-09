# AGENTS.md

Instructions for AI coding assistants and autonomous agents working in this repository. These apply regardless of which tool or model is being used.

## Do not modify this file

Pull requests that add, remove, or change `AGENTS.md` will not be reviewed or merged. If you believe a change is needed, open an issue for discussion instead.

## Tool-specific configuration

This file is the only assistant instruction file tracked in the repository. Tool-specific context or configuration files (for example `CLAUDE.md`, `GEMINI.md`, `.cursorrules`, or `.claude/` and similar directories) must not be committed. They are not listed in `.gitignore`; use your tool's user-level configuration, local ignore rules, or memory to keep such context outside the repository.

## Read and follow the contributor guidelines

Before making any change, read `CONTRIBUTING.rst` in the repository root and follow it in full. In particular:

- **Requirements for Code Submission**: tests for bug fixes and new features, documentation for new features, and passing pre-commit checks.
- **Generative AI Policy**: the human operating you is accountable for every change, must be able to explain every changed line, and must complete the AI Disclosure section in `.github/PULL_REQUEST_TEMPLATE.md`. Autonomous agents are not contributors. Do not open, comment on, or respond to review on pull requests or issues on a human's behalf.
- **Coding Style Guide**, **API Style Guide**, and **Docstrings**: match the existing conventions in the surrounding code.

Where this file and `CONTRIBUTING.rst` differ, `CONTRIBUTING.rst` wins.

## Code comments

- Do not add comments that describe what changed, why an edit was made, or which task or session produced it. Comments are not a changelog.
- Do not add comments restating what the code plainly does.
- Comment only when the logic itself would be non-obvious to an experienced developer reading it for the first time.
- Do not add, remove, or rewrite existing comments or docstrings outside the scope of the change you were asked to make.

## Commit messages

- Do not add AI attribution of any kind: no `Co-authored-by` lines for AI tools, no "Generated with", no tool names, no model names, and no links to session URLs.
- Write a short, descriptive subject line in the style of the existing history (`git log --oneline`).
- AI use is disclosed in the pull request, not in the commit log.

## Scope

- Make only the change that was requested. Do not refactor, reformat, or "clean up" unrelated code.
- Do not map, summarize, or document the codebase in this file or in new top-level files. Context about the code lives in the code, its docstrings, and the documentation under `doc/`; gather it as needed for the task at hand.
- Start with the tests closest to your change, then widen the pytest run to the affected subpackage, or the whole unit suite, when a change touches shared machinery. Tests that need sample data skip when no data directory is configured; that is expected. Do not attempt to run the nose-based answer tests unless explicitly requested.
