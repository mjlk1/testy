ALL=bin/bullet

.PHONY: all
all: $(ALL)

build/%.cpp.o: %.cpp
	@mkdir -p $(@D)
	g++ -mavx -c -Ofast -std=c++11 $< -o $@

SOURCES=$(shell find -wholename "./source/*.cpp")
OBJECTS=$(SOURCES:%=build/%.o)

bin/bullet: $(OBJECTS)
	@mkdir -p $(@D)
	g++ -s -flto $^ -o $@

%.out: %.in bin/bullet
	@mkdir -p $(@D)
	bin/bullet < $< > $@

TESTIN=$(shell find -wholename "./data/*.in")
TESTOUT=$(TESTIN:%.in=%.out)

.PHONY: test
test: $(TESTOUT)

.PHONY: clean
clean:
	rm -rRf build
	rm -f $(ALL)
	find data/ -name "*.out" -exec rm -f {} \;

