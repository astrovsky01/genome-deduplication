SHELL := /bin/sh

SYNOPSIS_SCRIPT := $(CURDIR)/code/synopsis.sh
DEFAULT_PREFIX := $(shell \
	if [ -d /opt/homebrew/bin ] && [ -w /opt/homebrew/bin ]; then \
		printf '%s' /opt/homebrew; \
	elif [ -d /usr/local/bin ] && [ -w /usr/local/bin ]; then \
		printf '%s' /usr/local; \
	else \
		printf '%s' "$(HOME)/.local"; \
	fi)
PREFIX ?= $(DEFAULT_PREFIX)
BIN_DIR := $(PREFIX)/bin
SYNOPSIS_BINARY := $(BIN_DIR)/synopsis

.DEFAULT_GOAL := all

.PHONY: all clean help FORCE

all: $(SYNOPSIS_BINARY)

$(SYNOPSIS_BINARY): FORCE $(SYNOPSIS_SCRIPT)
	@mkdir -p "$(BIN_DIR)"
	@rm -f "$@"
	@printf '%s\n' \
		'#!/bin/sh' \
		'if command -v zsh >/dev/null 2>&1; then' \
		'    exec zsh "$(SYNOPSIS_SCRIPT)" "$$@"' \
		'fi' \
		'exec bash "$(SYNOPSIS_SCRIPT)" "$$@"' > "$@"
	@chmod +x "$@"
	@printf 'Installed synopsis at %s\n' "$@"
	@case ":$$PATH:" in *":$(BIN_DIR):"*) ;; *) printf 'Add this to your shell profile: export PATH="%s:$$PATH"\n' "$(BIN_DIR)" ;; esac

FORCE:

clean:
	@rm -f "$(SYNOPSIS_BINARY)"

help:
	@printf '%s\n' \
		'Usage:' \
		'  make                       Install synopsis in a standard bin directory' \
		'  make PREFIX=/usr/local     Install to a specific prefix' \
		'  make clean                 Remove the installed executable'
