CC = g++
CCFLAGS += "-std=gnu++11"
CCFLAGS += "-fext-numeric-literals"
LDFLAGS += "-lquadmath"

stieltjes: main.o printer.o
	$(CC) $(CCFLAGS) main.o printer.o -o stieltjes -lquadmath

%: %.C
	$(CC) $(CCFLAGS) $< -o $@

%.o: %.C
	$(CC) $(CCFLAGS) -c $< -o $@

main.o: main.C OptionParser/Option.h OptionParser/OptionParser.h \
        printer.h stieltjes.h ql.h interpolant.h stieltjes_impl.h \
        input/InputReader.h
	$(CC) $(CCFLAGS) -c $< -o $@

printer.o: printer.C printer.h
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm -f *.o stieltjes
