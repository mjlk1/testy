ALL=bin/bullet

.PHONY: all
all: $(ALL)

build/%.cpp.o: %.cpp
	@mkdir -p $(@D)
	g++ -c $< -o $@

SOURCES=$(shell find -wholename "./source/*.cpp")
OBJECTS=$(SOURCES:%=build/%.o)

bin/bullet: $(OBJECTS)
	@mkdir -p $(@D)
	g++ $^ -o $@

.PHONY: clean
clean:
	rm -rRf build
	rm -f $(ALL)
