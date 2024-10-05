FROM node:8-alpine

RUN apk update && apk add --no-cache git python python-dev make gcc g++

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

ADD . /usr/src/app

RUN npm install -g serve
RUN npm install
RUN npm run build

CMD [ "serve", "-s", "build", "-l", "80" ]
