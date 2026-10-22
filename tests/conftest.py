# Copyright 2026 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Shared fixtures. The panel factory itself lives in ``tests/_panels.py``.

``LONGBET_TEST_VARIANT='{"linear_trt": true}'`` makes those keys defaults of
every ``LongBetConfig`` the tests build (explicit arguments still win), so a
test module can be rerun with the encourage-v2 options on.
"""

import ctypes
import gc
import json
import os

import pytest

if os.environ.get("LONGBET_TEST_VARIANT"):
    from longbet import _config as _cfg_mod

    _VARIANT = json.loads(os.environ["LONGBET_TEST_VARIANT"])
    _orig_init = _cfg_mod.LongBetConfig.__init__

    def _variant_init(self, *args, **kwargs):
        for key, value in _VARIANT.items():
            kwargs.setdefault(key, value)
        _orig_init(self, *args, **kwargs)

    _cfg_mod.LongBetConfig.__init__ = _variant_init

from _panels import make_staggered_panel


@pytest.fixture
def staggered_panel():
    return make_staggered_panel()


@pytest.fixture
def unbalanced_panel():
    return make_staggered_panel(seed=123, missing_rate=0.12)


@pytest.fixture
def binary_panel():
    return make_staggered_panel(seed=7, binary=True)


_MODULES_DONE = 0
_RSS_LIMIT_MODULE_MB = 2800.0
_RSS_LIMIT_TEST_MB = 3500.0
_LAST_CLEAR_RSS_MB = 0.0


def _current_rss_mb() -> float:
    """Return current process resident set size (VmRSS) in MiB on Linux."""
    try:
        with open("/proc/self/status", encoding="utf-8") as status_file:
            for line in status_file:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return 0.0


def _clear_jax_and_gc() -> None:
    global _LAST_CLEAR_RSS_MB
    import jax

    jax.clear_caches()
    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except OSError:
        pass
    _LAST_CLEAR_RSS_MB = _current_rss_mb()


@pytest.fixture(autouse=True)
def _guard_test_rss():
    yield
    if _current_rss_mb() >= max(_RSS_LIMIT_TEST_MB, _LAST_CLEAR_RSS_MB + 1200.0):
        _clear_jax_and_gc()


@pytest.fixture(autouse=True, scope="module")
def _clear_jax_compilation_cache():
    global _MODULES_DONE
    yield
    _MODULES_DONE += 1
    if _MODULES_DONE % 10 == 0 or _current_rss_mb() >= max(
        _RSS_LIMIT_MODULE_MB, _LAST_CLEAR_RSS_MB + 700.0
    ):
        _clear_jax_and_gc()

