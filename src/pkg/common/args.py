def remove_none(d: dict) -> dict:
    return {key: value for key, value in d.items() if value is not None}
