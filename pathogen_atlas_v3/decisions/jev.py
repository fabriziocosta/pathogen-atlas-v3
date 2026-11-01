"""OpenRouter typed Decisions API adapter for Jev System 1."""

from __future__ import annotations

import os
import time
from threading import local
from typing import Any, Protocol

import requests

DEFAULT_JEV_MODEL = "typesafe/jev-1.14"
DEFAULT_DECISIONS_URL = "https://openrouter.ai/api/alpha/decisions"


class DecisionClient(Protocol):
    model: str
    provider: str
    def decide_many(self, state: dict[str, Any], *, questions: dict[str, dict[str, Any]]) -> dict[str, Any]: ...


class JevRequestError(RuntimeError):
    def __init__(self, message: str, *, attempts: int):
        super().__init__(message)
        self.attempts = attempts


class OpenRouterJevClient:
    def __init__(self, *, api_key: str | None = None, model: str = DEFAULT_JEV_MODEL, endpoint: str = DEFAULT_DECISIONS_URL, timeout: float = 90, retries: int = 3, backoff_seconds: float = 1, session: requests.Session | None = None):
        self.api_key = (api_key or os.environ.get("OPENROUTER_API_KEY", "")).strip()
        if not self.api_key:
            raise ValueError("Set OPENROUTER_API_KEY before running Jev decisions.")
        if retries < 0 or timeout <= 0 or backoff_seconds < 0:
            raise ValueError("Invalid retry, timeout, or backoff setting.")
        self.model, self.provider, self.endpoint = str(model), "openrouter_jev", str(endpoint)
        self.timeout, self.retries, self.backoff_seconds = float(timeout), int(retries), float(backoff_seconds)
        self._provided_session = session
        self._thread_state = local()

    def remaining_balance(self) -> tuple[float | None, str]:
        """Best-effort remaining credit/limit lookup; never returns the secret key."""
        headers = {"Authorization": f"Bearer {self.api_key}"}
        for url, field in (("https://openrouter.ai/api/v1/credits", "credits"), ("https://openrouter.ai/api/v1/key", "key")):
            try:
                response = requests.get(url, headers=headers, timeout=8)
                if not response.ok:
                    continue
                data = response.json().get("data", {})
                if field == "credits":
                    total, used = data.get("total_credits"), data.get("total_usage")
                    if isinstance(total, (int, float)) and isinstance(used, (int, float)):
                        return max(0.0, float(total) - float(used)), "account credits"
                else:
                    left = data.get("limit_remaining")
                    if isinstance(left, (int, float)):
                        return max(0.0, float(left)), "API key limit"
            except Exception:
                continue
        return None, "remaining balance unavailable"

    @property
    def last_attempts(self) -> int:
        return int(getattr(self._thread_state, "last_attempts", 1))

    def decide_many(self, state: dict[str, Any], *, questions: dict[str, dict[str, Any]]) -> dict[str, Any]:
        if not questions:
            raise ValueError("At least one typed question is required.")
        payload = {"model": self.model, "state": state, "questions": {key: {"type": q["type"], "instructions": q["instructions"], "criteria": q.get("criteria", {})} for key, q in questions.items()}}
        headers = {"Authorization": f"Bearer {self.api_key}", "Content-Type": "application/json"}
        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            session = self._provided_session or getattr(self._thread_state, "session", None)
            if session is None:
                session = requests.Session()
                self._thread_state.session = session
            try:
                response = session.post(self.endpoint, json=payload, headers=headers, timeout=self.timeout)
                response.raise_for_status()
                value = response.json()
                if not isinstance(value, dict):
                    raise ValueError("Jev response must be a JSON object.")
                self._thread_state.last_attempts = attempt + 1
                return value
            except Exception as exc:
                last_error = exc
                transient = isinstance(exc, (requests.Timeout, requests.ConnectionError))
                if isinstance(exc, requests.HTTPError):
                    status = exc.response.status_code if exc.response is not None else 0
                    transient = status == 429 or 500 <= status < 600
                if not transient or attempt == self.retries:
                    break
                time.sleep(self.backoff_seconds * (2 ** attempt))
        raise JevRequestError(f"Jev request failed: {last_error}", attempts=attempt + 1) from last_error


def extract_noul_probability(response: dict[str, Any], question_id: str) -> float:
    answers = response.get("answers")
    answer = answers.get(question_id) if isinstance(answers, dict) else None
    if not isinstance(answer, dict) or answer.get("type") != "noul":
        raise ValueError(f"Missing {question_id!r} noul answer.")
    value = answer.get("noul")
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not 0 <= float(value) <= 1:
        raise ValueError(f"Invalid noul probability for {question_id!r}: {value!r}")
    return float(value)


def pathogen_gate_question(domain: str, prompt: dict[str, Any]) -> dict[str, Any]:
    question = prompt["domains"][domain]
    return {"human_pathogen_evidence": {"type": "noul", "instructions": question["instructions"], "criteria": question["criteria"]}}
