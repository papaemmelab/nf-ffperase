
FROM debian:bullseye-slim

ENV BIN_DIR=/usr/local/bin
RUN mkdir -p $BIN_DIR
COPY ./bin/annotate_w_pileup $BIN_DIR/
RUN chmod +x $BIN_DIR/annotate_w_pileup
ENV PATH=$BIN_DIR:$PATH

ENTRYPOINT ["annotate_w_pileup"]