SOURCES = $(wildcard *.Rnw)
FINALS = $(patsubst %.Rnw, %.pdf, $(SOURCES))

all: $(FINALS)

%.pdf: %.Rnw
	Sweave.sh -ld $<
