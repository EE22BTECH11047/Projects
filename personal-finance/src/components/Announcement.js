import React, { Component } from 'react';

class Announcement extends Component {
  render() {
    return (
      <section className="space-xs text-center bg-gradient text-light">
        <div className="container">
          <div className="row">
            <div className="col">
              <i className={this.props.icon + ' mr-1'} />
              <span className="mr-2">{this.props.text}</span>
              <a
                href={this.props.linkAddress}
                className="text-white"
                target="_blank"
                rel="noopener noreferrer"
              >
                {this.props.linkText}
              </a>
            </div>
          </div>
        </div>
      </section>
    );
  }
}

export default Announcement;
