# -*- coding: utf-8 -*-
import json
import logging

import requests

from ot_simple_connector.jobs import Jobs


class Connector:

    auth_ep = '/api/auth/login'

    def __init__(self, host, port, user, password, ssl=False, loglevel="INFO"):
        self.ssl = ssl
        self.password = password
        self.user = user
        self.port = port
        self.host = host

        self._session = None
        self._base_url = None

        logger = logging.getLogger('ot_simple_connector')
        ch = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s %(levelname)-s PID=%(process)d %(module)s:%(lineno)d func=%(funcName)s - %(message)s'
        )
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        logger.setLevel(loglevel)
        self.logger = logger

    @property
    def jobs(self):
        return Jobs(self.session)

    @property
    def base_url(self):
        if self._base_url is None:
            protocol = 'https' if self.ssl else 'http'
            url = '%s://%s:%s' % (protocol, self.host, self.port)
            return url
        else:
            return self._base_url

    @property
    def session(self):
        if self._session is None:
            self.logger.debug("Trying to connect...")
            session = requests.sessions.Session()
            response = session.post(
                self.base_url + self.auth_ep,
                data=json.dumps({'username': self.user, 'password': self.password})
            )
            self.logger.debug("Response.status_code: %s." % response.status_code)
            if response.status_code == 200:
                self.logger.debug("Response.text: %s." % response.text)
                response_json = response.json()
                response_status = response_json.get('status')
                if response_status == 'success':
                    self._session = session
                    self._session.base_url = self.base_url
                    self._session.authorized = True
                    self._session.username = self.user
                else:
                    raise Exception("Auth response: %s" % response.text)
            else:
                raise Exception("Auth code: %s" % response.status_code)

        return self._session
