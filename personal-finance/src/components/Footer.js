import React, { Component } from 'react';

class Footer extends Component {
  render() {
    return (
      <div className="container">
        <hr />
        <div className="text-center">
          <p className="text-small">
            &copy; 2020{' '}
            . All rights reserved.
          </p>
          <p className="text-small">
            Your data is saved only in the browser. If you switch browsers or
            clear the cache, it will not be persisted.
          </p>
        </div>
        <br />
      </div>
    );
  }
}

export default Footer;
