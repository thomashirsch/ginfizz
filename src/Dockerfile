
FROM miykael/nipype_level5


RUN echo "--- START   ----------"

# on monte un volume /rstp_data pour récupérer détarer le dataset et un autre pour mettre le code

VOLUME /rstp_data

VOLUME /rstp_code

ENV DATAROOT=/rstp_data

ENV CODEROOT=/rstp_code

COPY ./ginfizz_wrapper.sh  /rstp_code/ginfizz_wrapper.sh
RUN   /bin/chmod 777  /rstp_code/ginfizz_wrapper.sh

# on va recuperer le pipeline de preprocess en python

COPY ./ginfizz_preprocess.py  /rstp_code/ginfizz_preprocess.py


ENTRYPOINT ["/bin/bash"]



