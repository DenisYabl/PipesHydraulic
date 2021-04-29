# -*- coding: utf-8 -*-
import logging

from ot_simple_connector.dataset import Dataset


class Job:

    mj_ep = '/api/makejob'
    gr_ep = '/api/getresult'
    cj_ep = '/api/checkjob'

    logger = logging.getLogger('ot_simple_connector')

    def __init__(self, session, query_text, cache_ttl, tws, twf, sid):
        self.session = session
        self.twf = twf
        self.tws = tws
        self.cache_ttl = cache_ttl
        self.query_text = query_text
        self.sid = sid
        self.msg = None
        self.cid = None
        self._dataset = None

    @property
    def payload(self):
        _payload = {
            'original_otl': self.query_text,
            'cache_ttl': self.cache_ttl,
            'tws': self.tws,
            'twf': self.twf,
            'username': self.session.username,
            'sid': self.sid,
            'cid': self.cid
        }
        return _payload

    def create(self):
        response = self.session.post(self.session.base_url + self.mj_ep, data=self.payload)
        if response.status_code == 200:
            self.logger.debug('Response.text: %s.' % response.text)
            response_json = response.json()
            response_status = response_json.get('status')
            if response_status != 'success':
                raise Exception('Create error: %s' % response_status)
        else:
            raise Exception('Create http error: %s' % response.status_code)

    @property
    def status(self):
        response = self.session.get(self.session.base_url + self.cj_ep, params=self.payload)
        if response.status_code == 200:
            self.logger.debug('Response.text: %s.' % response.text)
            response_json = response.json()
            response_status = response_json.get('status')
            status = response_status
            self.msg = response_json.get('error')
            self.cid = response_json.get('cid')
        else:
            status = 'http_error'
            self.msg = 'http code: %s' % response.status_code
        return status

    @property
    def dataset(self):
        if self._dataset is None:
            if self.cid is None:
                raise Exception('Job with status %s has no cache id' % self.status)
            else:
                response = self.session.get(self.session.base_url + self.gr_ep, params=self.payload)
                self.logger.debug('Response code: %s.' % response)
                if response.status_code == 200:
                    self.logger.debug('Response.text: %s.' % response.text)
                    response_json = response.json()
                    response_status = response_json.get('status')
                    if response_status == 'success':
                        data_urls = response_json.get('data_urls')
                        urls = []
                        schema_url = None
                        for url in data_urls:
                            url = self.session.base_url + '/' + url
                            if '_SCHEMA' in url:
                                schema_url = url
                            else:
                                urls.append(url)
                        self.logger.debug('Schema url: %s.' % schema_url)
                        self.logger.debug('Data urls: %s.' % urls)
                        self._dataset = Dataset(self.session, schema_url, urls)
                    else:
                        raise Exception('Job with status %s has no dataset' % self.status)
                else:
                    raise Exception('Job with status %s has no dataset because of http error: %s' % (
                        self.status, response.status_code
                    ))

        return self._dataset
