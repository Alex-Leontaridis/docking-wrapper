from importlib import import_module

def create_component(kind: str, name: str, cfg: dict):
    """
    Dynamically import and instantiate a component.
    Example: kind='engines', name='vina' → engines.vina.VinaEngine
    """
    module = import_module(f"{kind}.{name}")
    cls = getattr(module, name.capitalize() + kind.capitalize())
    return cls(**cfg) 