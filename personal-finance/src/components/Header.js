import React, { Component } from 'react';
import NumberFormat from 'react-number-format';
import ReactTooltip from 'react-tooltip';

import CONFIG from '../config';

class Header extends Component {
  constructor(props) {
    super(props);

    this.state = {
      accrued_sum: 0,
      monthly_income: 0,
      monthly_expenses: 0,
      monthly_total: 0,
    };

    this.interval = null;
  }

  computeData() {
    let monthly_income = 0;
    let monthly_expenses = 0;

    this.props.items.map(function (item) {
      if (item.type === CONFIG.ITEM_TYPE.ASSET) {
        monthly_income += (item.amount * item.apr) / 100 / 12;
      } else if (item.type === CONFIG.ITEM_TYPE.LIABILITY) {
        monthly_expenses += (item.amount * item.apr) / 100 / 12;
      } else if (item.type === CONFIG.ITEM_TYPE.INCOME) {
        monthly_income += item.amount;
      } else if (item.type === CONFIG.ITEM_TYPE.EXPENSE) {
        monthly_expenses += item.amount;
      }

      return item;
    });

    this.setState({
      accrued_sum: 0,
      monthly_income: monthly_income,
      monthly_expenses: monthly_expenses,
      monthly_total: monthly_income - monthly_expenses,
    });
  }

  componentDidMount() {
    const _this = this;
    let accrued_sum = 0;

    this.computeData();

    this.interval = setInterval(function () {
      _this.props.items.map(function (item) {
        let accrued = 0;

        if (item.type === CONFIG.ITEM_TYPE.ASSET) {
          accrued +=
            ((item.amount * item.apr) / 100 / CONFIG.MILISECONDS_IN_YEAR) *
            CONFIG.REFRESH_INTERVAL;
        } else if (item.type === CONFIG.ITEM_TYPE.LIABILITY) {
          accrued -=
            ((item.amount * item.apr) / 100 / CONFIG.MILISECONDS_IN_YEAR) *
            CONFIG.REFRESH_INTERVAL;
        } else if (item.type === CONFIG.ITEM_TYPE.INCOME) {
          accrued +=
            ((item.amount * 12) / CONFIG.MILISECONDS_IN_YEAR) *
            CONFIG.REFRESH_INTERVAL;
        } else if (item.type === CONFIG.ITEM_TYPE.EXPENSE) {
          accrued -=
            ((item.amount * 12) / CONFIG.MILISECONDS_IN_YEAR) *
            CONFIG.REFRESH_INTERVAL;
        }

        accrued_sum += accrued;

        return item;
      });

      _this.setState({
        accrued_sum: accrued_sum,
      });
    }, CONFIG.REFRESH_INTERVAL);
  }

  componentWillUnmount() {
    if (this.interval) {
      clearInterval(this.interval);
    }
  }

  render() {
    return (
      <div className="jumbotron jumbotron-fluid">
        <div className="container">
          <div className="row">
            <ReactTooltip effect="solid" />

            <div className="col-5">
              <h1 className="display-4">Money Streams</h1>
              <p className="lead">Keep track of income vs expenses and view diagrams &amp; charts to help you get a grasp on your financial situation.</p>
              <button
                type="button"
                className={`btn btn-lg btn-block ${
                  this.state.accrued_sum >= 0 ? 'btn-primary' : 'btn-danger'
                }`}
                data-tip
                data-for="streaming_total"
              >
                Streaming&nbsp;&nbsp;
                <h5 style={{ color: 'rgb(177, 190, 212)' }}>
                  <NumberFormat
                    value={this.state.accrued_sum}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={6}
                  />
                </h5>
                <ReactTooltip id="streaming_total" effect="solid">
                  Hourly:{' '}
                  <NumberFormat
                    value={this.state.monthly_total / 30 / 24}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={2}
                  />
                  <br />
                  Daily:{' '}
                  <NumberFormat
                    value={this.state.monthly_total / 30}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={1}
                  />
                  <br />
                  <strong>
                    Monthly:{' '}
                    <NumberFormat
                      value={this.state.monthly_total}
                      displayType={'text'}
                      thousandSeparator={true}
                      prefix={'$'}
                      decimalScale={0}
                    />
                  </strong>
                  <br />
                  Yearly:{' '}
                  <NumberFormat
                    value={this.state.monthly_total * 12}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={0}
                  />
                </ReactTooltip>
              </button>
              <br />
              <br />
              <h5>
                <strong>Total Income:</strong>&nbsp;
                <span className="badge badge-primary" data-tip="Hourly">
                  <NumberFormat
                    value={this.state.monthly_income / 30 / 24}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={2}
                  />
                </span>{' '}
                /&nbsp;
                <span className="badge badge-primary" data-tip="Daily">
                  <NumberFormat
                    value={this.state.monthly_income / 30}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={1}
                  />
                </span>{' '}
                /&nbsp;
                <span className="badge badge-primary" data-tip="Monthly">
                  <strong>
                    <NumberFormat
                      value={this.state.monthly_income}
                      displayType={'text'}
                      thousandSeparator={true}
                      prefix={'$'}
                      decimalScale={0}
                    />
                  </strong>
                </span>{' '}
                /&nbsp;
                <span className="badge badge-primary" data-tip="Yearly">
                  <NumberFormat
                    value={this.state.monthly_income * 12}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={0}
                  />
                </span>
              </h5>
              <p>
                <strong>Expenses:</strong>&nbsp;
                <span className="badge badge-danger" data-tip="Hourly">
                  <NumberFormat
                    value={this.state.monthly_expenses / 30 / 24}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={2}
                  />
                </span>{' '}
                /&nbsp;
                <span className="badge badge-danger" data-tip="Daily">
                  <NumberFormat
                    value={this.state.monthly_expenses / 30}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={1}
                  />
                </span>{' '}
                /&nbsp;
                <span className="badge badge-danger" data-tip="Monthly">
                  <strong>
                    <NumberFormat
                      value={this.state.monthly_expenses}
                      displayType={'text'}
                      thousandSeparator={true}
                      prefix={'$'}
                      decimalScale={0}
                    />
                  </strong>
                </span>{' '}
                /&nbsp;
                <span className="badge badge-danger" data-tip="Yearly">
                  <NumberFormat
                    value={this.state.monthly_expenses * 12}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={0}
                  />
                </span>
              </p>
            </div>
            <div className="col-sm">{this.props.children}</div>
          </div>
          <div className="row">
            <div className="col">
              <center>
                <button
                  className="btn btn-sm btn-secondary"
                  onClick={this.props.sampleHandler}
                  data-tip="Temporarily load sample data. Refresh page to reset."
                >
                  Load sample data
                </button>
              </center>
            </div>
          </div>
        </div>
      </div>
    );
  }
}

export default Header;
