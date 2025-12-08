from datetime import datetime
from typing import Any, Dict, Iterable, Type, Union
import requests
import json
from biggr import models

API_URL = "http://localhost:80/api/v3/"
OBJECTS_API_URL = f"{API_URL}objects/"
IDENTIFIERS_API_URL = f"{API_URL}identifiers/"

MODEL_NAMES = {x.__name__: x for x in models.Base.__subclasses__()}


def _request(api_url: str, data: Dict[str, Any]):
    r = requests.post(api_url, json=data)
    if r.status_code != 200:
        print(f"Status code: {r.status_code}")
        print(data)
        print(r.response)
        return None
    return r.json()


def get_raw(obj_type: Union[str, Type[models.Base]], obj_id: Union[str, int]):
    if not isinstance(obj_type, str):
        obj_type = obj_type.__name__
    data = {"type": obj_type, "id": obj_id}
    return _request(OBJECTS_API_URL, data)


def convert_result_to_models(o):
    if isinstance(o, str):
        return o
    if isinstance(o, list):
        return [convert_result_to_models(x) for x in o]
    if isinstance(o, dict):
        cls_type = dict
        if "_type" in o:
            if o["_type"] == "datetime":
                return datetime.fromisoformat(o["iso"])
            cls_type = MODEL_NAMES.get(o["_type"])
            if cls_type is None:
                raise ValueError()
        entries = {k: convert_result_to_models(v) for k, v in o.items() if k != "_type"}
        if cls_type is dict:
            return entries
        return cls_type.from_dict(entries)
    return o


def get(obj_type: Union[str, Type[models.Base]], obj_id: Union[str, int]):
    print(f"GET: {obj_type}: {obj_id}")
    result = get_raw(obj_type, obj_id)
    if result is None:
        return None
    if "object" in result:
        return convert_result_to_models(result["object"])
    else:
        return convert_result_to_models(result["objects"])


def get_metabolites_by_identifiers(
    identifiers: Union[str, Iterable], model_bigg_id: Optional[str] = None
):
    if isinstance(identifiers, str):
        identifiers = [identifiers]
    for identifier in identifiers:
        if not ":" in identifier:
            raise ValueError("Identifiers should be supplied as '<namespace>:<id>'.")
    query = {
        "type": "metabolite",
        "identifiers": identifiers,
        "model_bigg_id": model_bigg_id,
    }
    result = _request(IDENTIFIERS_API_URL, query)
    return convert_result_to_models(result)
