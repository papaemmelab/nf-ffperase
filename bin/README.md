# annotate_w_pileup

## Compile binary

The get the static binary, the script needs to be compiled using nim.

1. Install [`nim`](<https://nim-lang.org/install_unix.html>). Example, for MacOS: `brew install nim`

2. Get the `hts-nim_static_builder` binary.

    You can either download it from [here](https://github.com/brentp/hts-nim?tab=readme-ov-file#static-builds), or compile it (needs `docker`):

    ```bash
    git clone https://github.com/brentp/hts-nim 
    cd hts-nim/docker
    nim c -r "hts_nim_static_builder.nim" -d:release
    ```

3. Use the static builder to compile `annotate_w_pileup.nim`:

    ```bash
    ${path-to-binary}/hts_nim_static_builder -s scripts/annotate_w_pileup.nim --deps "hts@>=0.2.13" --deps "https://github.com/brentp/hileup.git" --deps "progress"
    ```

### For Apple M1/M2 users

The static builder does not work on Apple M1/M2 chips. You can compile the script using the following command: