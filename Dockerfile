FROM python:3.6.5

# Fix APT sources to use Debian Archive
RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list && \
    echo "deb http://archive.debian.org/debian-security stretch/updates main" >> /etc/apt/sources.list && \
    echo "Acquire::Check-Valid-Until false;" > /etc/apt/apt.conf.d/99no-check-valid-until && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    zlib1g-dev \
    libssl-dev \
    libffi-dev \
    libbz2-dev \
    libsqlite3-dev \
    curl \
    ca-certificates \
    git \
    autoconf \
    automake \
    && rm -rf /var/lib/apt/lists/*

# Set up locales
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# Install Picard 2.9.0
RUN apt-get update && apt-get install -y openjdk-8-jdk && \
    mkdir -p /downloads && \
    wget -qO /downloads/picard.jar "https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar"
ENV PICARD=/downloads/picard.jar

# Install HTSlib using libdeflate
RUN \
    cd /downloads && \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git --branch v1.2 && \
    cd libdeflate && \
    make -j$(nproc) CFLAGS='-fPIC -O3' libdeflate.a && \
    cp libdeflate.a /usr/local/lib && \
    cp libdeflate.h /usr/include && \
    cd /downloads && rm -rf libdeflate && \
    git clone https://github.com/samtools/htslib --branch 1.9 && \
    cd htslib && \
    autoheader && autoconf && \
    ./configure --enable-plugins --with-libdeflate && \
    make -j$(nproc) install && \
    cd /downloads && rm -rf htslib

# Install BEDOPS tools for vcf2bed
RUN \
    cd /downloads && \
    wget -q "https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2" && \
    tar jxf bedops_linux_x86_64-v2.4.40.tar.bz2 && \
    cp bin/* /usr/local/bin && \
    rm -rf bedops_linux_x86_64-v2.4.40.tar.bz2 bin

# Install split_bed_by_index
RUN wget -qO /usr/bin/split_bed_by_index \
    "https://github.com/papaemmelab/split_bed_by_index/releases/download/0.2.0b/split_bed_by_index" && \
    chmod +x /usr/bin/split_bed_by_index

# Install required Python packages
COPY requirements.txt /tmp/
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# Set up binary directory
ENV BIN_DIR=/usr/local/bin
RUN mkdir -p $BIN_DIR
COPY ./bin/annotate_w_pileup $BIN_DIR/
RUN chmod +x $BIN_DIR/annotate_w_pileup
ENV PATH=$BIN_DIR:$PATH

ENTRYPOINT ["annotate_w_pileup"]
