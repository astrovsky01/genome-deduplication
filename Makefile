SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2
LDFLAGS ?=
LDLIBS ?= -lz

SOURCE_DIR := $(CURDIR)/code
DEDUP_SOURCE := $(SOURCE_DIR)/dedup6.cpp
DEDUP_BINARY := $(SOURCE_DIR)/dedup6
DEFAULT_PREFIX := $(HOME)/.local
PREFIX ?= $(DEFAULT_PREFIX)
BIN_DIR := $(PREFIX)/bin
RUNTIME_DIR := $(PREFIX)/lib/synopsis
SCRIPT_SOURCES := $(wildcard $(SOURCE_DIR)/*.sh)
SYNOPSIS_BINARY := $(BIN_DIR)/synopsis

.DEFAULT_GOAL := all

.PHONY: all clean help install FORCE

all: $(DEDUP_BINARY) $(SYNOPSIS_BINARY) $(RUNTIME_DIR)

$(DEDUP_BINARY): $(DEDUP_SOURCE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 -o "$@" "$<" $(LDFLAGS) $(LDLIBS)
	@chmod +x "$@"
	@printf 'Built dedup6 at %s\n' "$@"

$(SYNOPSIS_BINARY): FORCE
	@mkdir -p "$(BIN_DIR)"
	@rm -f "$@"
	@printf '%s\n' \
		'#!/bin/sh' \
		'launcher_dir=$$(CDPATH= cd -- "$$(dirname -- "$$0")" && pwd)' \
		'if command -v zsh >/dev/null 2>&1; then' \
		'    exec zsh "$$launcher_dir/../lib/synopsis/synopsis.sh" "$$@"' \
		'fi' \
		'exec bash "$$launcher_dir/../lib/synopsis/synopsis.sh" "$$@"' > "$@"
	@chmod +x "$@"
	@printf 'Installed synopsis at %s\n' "$@"
	@case ":$$PATH:" in *":$(BIN_DIR):"*) ;; *) printf 'Add this to your shell profile: export PATH="%s:$$PATH"\n' "$(BIN_DIR)" ;; esac

$(RUNTIME_DIR): $(DEDUP_BINARY) $(SCRIPT_SOURCES)
	@mkdir -p "$@"
	@cp $(SOURCE_DIR)/*.sh "$@/"
	@cp "$(DEDUP_BINARY)" "$@/dedup6"
	@chmod +x "$@"/*.sh "$@/dedup6"

install: all $(RUNTIME_DIR)

FORCE:

clean:
	@rm -f "$(SYNOPSIS_BINARY)" "$(DEDUP_BINARY)"

help:
	@printf '%s\n' \
		'Usage:' \
		'  make                       Build dedup6 and install synopsis' \
		'  make PREFIX=/usr/local     Install to a specific prefix' \
		'  make install PREFIX=/usr/local' \
		'  make clean                 Remove the installed executable'
