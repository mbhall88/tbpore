PROJECT := "tbpore"
OPEN := if os() == "macos" { "open" } else { "xdg-open" }
VERSION := `poetry version | grep -oE '\d+\.\d+\.\d+'`
BOLD := `tput bold`
NORM := `tput sgr0`

# format code with black and isort
fmt:
    poetry run black .
    poetry run isort .

# check format of code with black and isort
check-fmt:
    poetry run black --check .
    poetry run isort --check .

# lint code with flake8
lint:
    poetry run flake8 .

# install latest version with poetry
install:
    poetry install

# run all tests
test:
    poetry run pytest tests/

# run tests with coverage report
coverage:
    poetry run pytest --cov-report term --cov-report html --cov={{ PROJECT }} --cov-branch tests/
    {{ OPEN }} htmlcov/index.html

# run tests on the CI
test-ci:
    poetry run pytest --cov={{ PROJECT }} --cov-report=xml --cov-branch tests/

# check formatting, linting, and tests
check: check-fmt lint test

# prints out the commands to run to tag the release and push it
tag:
    @echo "Run {{ BOLD }}git tag -a {{ VERSION }} -m <message>{{ NORM }} to tag the release"
    @echo "Then run {{ BOLD }}git push origin {{ VERSION }}{{ NORM }} to push the tag"
