from typing import Any, Dict

def kwarg_mgr(kwargs: dict, key: str, default: Any) -> Any:
    return kwargs[key] if key in kwargs else default

def sns_color_mgr(
    kwargs: dict,
    default: Dict[str, Any] = {}
) -> Any:
    default_color = kwarg_mgr(default, "color", None)
    default_palette = kwarg_mgr(default, "palette", None)
    
    if "color" in kwargs:
        color = kwargs["color"]
        palette = None
    
    elif "palette" in kwargs:
        color = None
        palette = kwargs["palette"]
    
    else:
        color = default_color
        palette = default_palette
    
    return {"color": color, "palette": palette}