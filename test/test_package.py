"""Package-level smoke tests: public API and back-compat re-exports."""

import pyquiver


def test_public_api_exports():
    for name in pyquiver.__all__:
        assert hasattr(pyquiver, name), name


def test_version_string():
    v = pyquiver.__version__
    assert isinstance(v, str)
    # either a real version (installed) or the documented source-checkout fallback
    assert "." in v or v == "0+unknown"


def test_orca_backcompat_reexport():
    # pyquiver.orca.parse_orca_output must remain importable
    from pyquiver.orca import parse_orca_output, parse
    assert callable(parse_orca_output)
    assert callable(parse)


def test_settings_backcompat():
    from pyquiver import settings
    assert settings.DEBUG == 0
