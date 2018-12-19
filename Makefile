all: ssss-split ssss-combine ssss-test ssss.1 ssss.1.html

ifeq ($(shell uname),Darwin)
OSX_SDK:=$(shell xcrun --sdk macosx --show-sdk-path)
CFLAGS+=\
	-I/opt/manual/include -L/opt/manual/lib \
	-DNOMLOCK \
	-D__APPLE__ \
	-isysroot $(OSX_SDK)
endif

TEST_CFLAGS= \
	     -O0 \
	     -ggdb \
	     -fno-omit-frame-pointer \
	     -fsanitize=address,undefined,unsigned-integer-overflow \
	     -fstack-protector-strong

ssss-split: ssss.c
	$(CC) $(CFLAGS) -W -Wall -Werror -O2 -lgmp -lcrypto -o ssss-split ssss.c
	strip ssss-split

ssss-combine: ssss-split
	ln -f ssss-split ssss-combine

ssss-test: ssss-test.c ssss.c
	$(CC) $(CFLAGS) $(TEST_CFLAGS) -DTEST -W -Wall -Werror -lgmp -lcrypto -o ssss-test ssss.c ssss-test.c

ssss.1: ssss.manpage.xml
	xmltoman ssss.manpage.xml > ssss.1

ssss.1.html: ssss.manpage.xml
	xmlmantohtml ssss.manpage.xml > ssss.1.html

clean:
	rm -rf ssss-split ssss-combine ssss-test ssss.1 ssss.1.html
