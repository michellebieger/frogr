"""Test IO functions."""


def test_namelist():
    """Test namelist generation."""
    from atmopy.io import create_namelist

    expected = """&TEST
param1 = 12.4
param2 = 27
param3 = 'hello there'
param4 = .false.
/"""

    test_value = create_namelist(
        "test",
        {"param1": 12.4, "param2": 27, "param3": "hello there", "param4": False},
    )

    assert test_value == expected

    expected = """&TEST
param1 = 12.4
param2 = 27
param3 = 'hello there'
param4 = .true.
/"""

    test_value = create_namelist(
        "test",
        {"param1": 12.4, "param2": 27, "param3": "hello there", "param4": True},
    )

    assert test_value == expected

    expected = """&TEST2
/"""
    assert expected == create_namelist("test2")
