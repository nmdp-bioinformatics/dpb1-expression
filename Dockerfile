##
## Copyright (c) 2024 NMDP.
##
## This file is part of the NMDP DP Tool 
## (see https://github.com/nmdp-bioinformatics/dpb1-expression).
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.
##

FROM node:14-alpine as angular-builder

LABEL maintainer="NMDP Bioinformatics"

COPY webapp/package.json webapp/package-lock.json ./
RUN npm ci && mkdir /ng-app

WORKDIR /ng-app
COPY ./webapp/ .
RUN npm install -g @angular/cli && npm install
ARG CONFIGURATION=production
RUN npm run ng build -- --configuration=$CONFIGURATION

FROM python:3.11
#FROM dpb1-expression-tce-base-image:0.0.1-e18a8c78

COPY . /app
WORKDIR /app
RUN rm -r /app/webapp/*
RUN mkdir -p /app/logs/

RUN adduser -D -g dpb1tceexp -u 1000 dpb1tceexp && \
	chown -R dpb1tceexp:dpb1tceexp /app && \
  mkdir -p /var/log/ /var/lib/nginx/ /run/nginx && \
  chmod -R 777 /var/log/ /var/lib/nginx/ /run/nginx

# Upgrade pip

COPY --from=angular-builder ./ng-app/dist/dp_tool /usr/share/nginx/html
COPY webapp/nginx-conf/default.conf /etc/nginx/http.d/
COPY docker-entrypoint-unified-prod.sh /usr/local/bin/

EXPOSE 8080

ENTRYPOINT [ "/usr/local/bin/docker-entrypoint-unified-prod.sh" ]