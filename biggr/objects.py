from datetime import datetime
from typing import Any, Dict, Iterable, Optional, Type, Union
import requests
from biggr import models

API_URL = "https://biggr.org/api/v3/"
OBJECTS_API_URL = f"{API_URL}objects/"
IDENTIFIERS_API_URL = f"{API_URL}identifiers/"

MODEL_NAMES = {x.__name__: x for x in models.Base.__subclasses__()}


def _request(api_url: str, data: Dict[str, Any]) -> Optional[Any]:
    r = requests.post(api_url, json=data)
    if r.status_code != 200:
        print(f"Status code: {r.status_code}")
        print(data)
        return None
    return r.json()


def get_raw(
    obj_type: Union[str, Type[models.Base]], obj_id: Union[str, int]
) -> Optional[Any]:
    """Get the raw API request results as a python dictionary.
    
    Makes the BiGGr API request to obtain object(s) of type `obj_type`. Returns the JSON
    result interpreted as python dictionary.

    Parameters
    ----------
    obj_type: str or the class of the object to be retrieved
        For all possible values, please refer to the BiGGr API documentation at
        `biggr.org/data_access <https://biggr.org/data_access#data_objects_objects>`__
    obj_id: str or int
        If `obj_id` is a str, the ID is interpreted as a BiGG ID (this can not be used
        with all database entities). In the case that `obj_id` is of type int, the ID is
        interpreted as an internal ID, as used for defining relationships between
        database entities.
    """
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
    """Get an entity from the BiGGr database and return it as a python object.
    
    Makes the BiGGr API request to obtain object(s) of type `obj_type`. Returns the JSON
    result interpreted as an instance of the applicable class, as available in the
    `biggr.models` module. Some relations are loaded by default, whilst others are
    loaded automatically when accessed.

    Parameters
    ----------
    obj_type: str or the class of the object to be retrieved
        For all possible values, please refer to the BiGGr API documentation at
        `biggr.org/data_access <https://biggr.org/data_access#data_objects_objects>`__
    obj_id: str or int
        If `obj_id` is a str, the ID is interpreted as a BiGG ID (this can not be used
        with all database entities). In the case that `obj_id` is of type int, the ID is
        interpreted as an internal ID, as used for defining relationships between
        database entities.
    """
    # print(f"GET: {obj_type}: {obj_id}")
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
