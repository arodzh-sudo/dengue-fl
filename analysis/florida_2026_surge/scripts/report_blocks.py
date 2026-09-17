#!/usr/bin/env python3
"""The block vocabulary a build's content.py is written in."""


def block(kind, **fields):
    fields["kind"] = kind
    return fields


def read_config(path):
    """A small subset of YAML: key: value, inline lists, and one level of nesting."""
    config, section = {}, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.split("#", 1)[0].rstrip()
            if not line.strip():
                continue
            indented = line.startswith((" ", "\t"))
            key, _, value = line.strip().partition(":")
            key, value = key.strip(), value.strip()
            if indented and section is not None:
                config[section][key] = strip_quotes(value)
                continue
            if not value:
                section = key
                config[key] = {}
                continue
            section = None
            if value.startswith("[") and value.endswith("]"):
                config[key] = [strip_quotes(v) for v in value[1:-1].split(",") if v.strip()]
            else:
                config[key] = strip_quotes(value)
    return config


def strip_quotes(value):
    value = value.strip()
    if len(value) > 1 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value
