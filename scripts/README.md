# Shared Scripts

This folder contains shared utilities used by the HCV workflows. Each utility has its own folder and is run as `scripts/<name>/<name>.py`.

Use a script's built-in help for its inputs and outputs:

```bash
uv run python scripts/<script_name>/<script_name>.py --help
```

The scripts cover COMET profile post-processing, sequence audits, reference and consensus reporting, GenBank/reference maintenance, and filtering-support exports.
