# gatk-minimap2-jni
JNI code for minimap2

This project allows Java code to call Heng Li's minimap2 assembler.
We include binaries for OSX and Linux in a maven central artifact.

To build and install locally you'll need gmake, git, gcc, and Java 8. Execute
```
./gradlew install
```

This will work for testing, but will only produce a native library for your local system.

To use this JNI binding on an architecture other than OSX or Linux:
  Go into ```src/main/c```.
  Modify the Makefile to produce a library name appropriate to your system.
  Type ```make``` (you'll need gmake, git, and gcc).
  Move the library you built somewhere permanent on your machine.
  Use ```-DLIBMM2_PATH=<that permanent location>``` when you run GATK (or other Java program).
