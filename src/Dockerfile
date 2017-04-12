
FROM miykael/nipype_level2


RUN echo "--- START   ----------"

# 30 mars hackfest on monte une directory /rstp_data pour récupérer détarer le dataset et un autre pour mettre le code

RUN mkdir /rstp_data

RUN mkdir /rstp_code

COPY ./ginfizz_wrapper.sh  /rstp_code/ginfizz_wrapper.sh
RUN   chmod 777  /rstp_code/ginfizz_wrapper.sh

# on va recuperer le pipeline de preprocess en python

COPY ./ginfizz_preprocess.py  /rstp_code/ginfizz_preprocess.py


ENTRYPOINT ["/bin/bash"]



