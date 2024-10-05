import React, { Component } from 'react';
import NumberFormat from 'react-number-format';
import ReactTooltip from 'react-tooltip';

import CONFIG from '../config';

import Pie from './Pie';
import Storage from './Storage';

class Expenses extends Component {
  constructor(props) {
    super(props);

    this.interval = null;
    this.state = {
      accrued_sum: 0,
      name: '',
      expense: '',
      show_add_form: false,
    };

    this.handleInputChange = this.handleInputChange.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
    this.showAddForm = this.showAddForm.bind(this);
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.type === 'checkbox' ? target.checked : target.value;
    const name = target.name;

    this.setState({
      [name]: value,
    });
  }

  handleSubmit(e) {
    e.preventDefault();

    Storage.addItem({
      id: Date.now(),
      type: CONFIG.ITEM_TYPE.EXPENSE,
      name: this.state.name,
      amount: this.state.expense,
      apr: 0,
      compounds: false,
    });

    this.setState({
      name: '',
      expense: '',
    });

    window.location.reload();
  }

  handleDelete(item) {
    Storage.removeItem(item.id);

    window.location.reload();
  }

  showAddForm() {
    this.setState({
      show_add_form: true,
    });
  }

  componentDidMount() {
    const _this = this;

    this.interval = setInterval(function () {
      let accrued_sum = 0;

      _this.props.items.map(function (item) {
        let accrued = item.accrued || 0;

        accrued +=
          ((item.amount * 12) / CONFIG.MILISECONDS_IN_YEAR) *
          CONFIG.REFRESH_INTERVAL;

        item.accrued = accrued;

        accrued_sum += accrued;

        return item;
      });

      _this.setState({
        accrued_sum: accrued_sum,
      });

      _this.forceUpdate();
    }, CONFIG.REFRESH_INTERVAL);
  }

  componentWillUnmount() {
    if (this.interval) {
      clearInterval(this.interval);
    }
  }

  getTotalEstimatedPrincipal() {
    let total_principal = 0;

    this.props.items.forEach(function (asset) {
      const estimated_principal =
        ((asset.amount * 12) / CONFIG.INCOME_ESTIMATE_APR) * 100;
      total_principal += estimated_principal;
    });

    return total_principal + this.state.accrued_sum;
  }

  getMonthlyDebit() {
    let monthly_debit = 0;

    this.props.items.forEach(function (asset) {
      monthly_debit += asset.amount;
    });

    return monthly_debit;
  }

  render() {
    const { name, expense } = this.state;
    const valid = name.length > 0 && expense.length > 0;

    return (
      <React.Fragment>
        <ReactTooltip id="expenses_debit_tooltip" effect="solid">
          Hourly:{' '}
          <NumberFormat
            value={this.getMonthlyDebit() / 30 / 24}
            displayType={'text'}
            thousandSeparator={true}
            prefix={'$'}
            decimalScale={2}
          />
          <br />
          Daily:{' '}
          <NumberFormat
            value={this.getMonthlyDebit() / 30}
            displayType={'text'}
            thousandSeparator={true}
            prefix={'$'}
            decimalScale={1}
          />
          <br />
          <strong>
            Monthly:{' '}
            <NumberFormat
              value={this.getMonthlyDebit()}
              displayType={'text'}
              thousandSeparator={true}
              prefix={'$'}
              decimalScale={0}
            />
          </strong>
          <br />
          Yearly:{' '}
          <NumberFormat
            value={this.getMonthlyDebit() * 12}
            displayType={'text'}
            thousandSeparator={true}
            prefix={'$'}
            decimalScale={0}
          />
        </ReactTooltip>

        <ReactTooltip id="expenses_chart_tooltip" effect="solid" place="left">
          <Pie
            items={this.props.items}
            title="Expenses"
            backgroundColor="white"
          />
        </ReactTooltip>

        <div className="row">
          <div className="col-sm">
            <h3>
              <i className="icon-shopping-cart" />
              &nbsp;
              <strong data-tip data-for="expenses_chart_tooltip">
                Expenses
              </strong>
            </h3>
            <span data-tip data-for="expenses_debit_tooltip">
              Debit:{' '}
              <NumberFormat
                value={this.state.accrued_sum}
                displayType={'text'}
                thousandSeparator={true}
                prefix={'$'}
                decimalScale={5}
              />
            </span>
          </div>
        </div>

        <hr />

        <div
          className="row"
          style={{ display: this.state.show_add_form ? 'none' : 'block' }}
        >
          <div className="col">
            <center>
              <button className="btn btn-danger" onClick={this.showAddForm}>
                <i className="icon-circle-with-plus" />
              </button>
            </center>
          </div>
        </div>

        <div
          className="row"
          style={{ display: this.state.show_add_form ? 'block' : 'none' }}
        >
          <div className="col-sm">
            <form onSubmit={this.handleSubmit}>
              <div className="form-row form-group">
                <div className="col-6">
                  <input
                    className="form-control form-control-lg"
                    type="text"
                    placeholder="Name"
                    name="name"
                    value={this.state.name}
                    onChange={this.handleInputChange}
                  />
                  <small>eg. Taxes</small>
                </div>
                <div className="col-6">
                  <input
                    className="form-control form-control-lg"
                    type="number"
                    placeholder="Monthly Expense"
                    name="expense"
                    value={this.state.expense}
                    onChange={this.handleInputChange}
                  />
                  <small>eg. $1000</small>
                </div>
              </div>
              <div className="form-row form-group">
                <div className="col">
                  <button
                    className="btn btn-block btn-danger btn-lg"
                    disabled={!valid}
                  >
                    Add expense
                  </button>
                </div>
              </div>
            </form>
          </div>
        </div>

        <hr />

        <table className="table table-hover align-items-center table-borderless">
          <thead>
            <tr>
              <th scope="col">
                <strong>Expense</strong>
              </th>
              <th scope="col">
                <strong>Amount</strong>
              </th>
              <th scope="col">&nbsp;</th>
            </tr>
          </thead>
          {this.props.items.map((item) => (
            <tbody
              key={item.id}
              style={{ borderTop: 'none' }}
              data-tip
              data-for={`item_tooltip_${item.id}`}
            >
              <ReactTooltip id={`item_tooltip_${item.id}`} effect="solid">
                Hourly:{' '}
                <NumberFormat
                  value={item.amount / 30 / 24}
                  displayType={'text'}
                  thousandSeparator={true}
                  prefix={'$'}
                  decimalScale={2}
                />
                <br />
                Daily:{' '}
                <NumberFormat
                  value={item.amount / 30}
                  displayType={'text'}
                  thousandSeparator={true}
                  prefix={'$'}
                  decimalScale={1}
                />
                <br />
                Monthly:{' '}
                <NumberFormat
                  value={item.amount}
                  displayType={'text'}
                  thousandSeparator={true}
                  prefix={'$'}
                  decimalScale={0}
                />
                <br />
                Yearly:{' '}
                <NumberFormat
                  value={item.amount * 12}
                  displayType={'text'}
                  thousandSeparator={true}
                  prefix={'$'}
                  decimalScale={0}
                />
              </ReactTooltip>
              <tr className="bg-white" key={item.id}>
                <th>
                  <strong>{item.name}</strong>
                </th>
                <td>
                  <NumberFormat
                    value={item.amount}
                    displayType={'text'}
                    thousandSeparator={true}
                    prefix={'$'}
                    decimalScale={0}
                  />
                </td>
                <td>
                  <button
                    className="btn btn-secondary btn-sm"
                    onClick={this.handleDelete.bind(this, item)}
                  >
                    <i className="icon-circle-with-cross" />
                  </button>
                </td>
              </tr>
              <tr className="table-divider"></tr>
            </tbody>
          ))}
        </table>

        {this.props.items.length === 0 && <center>No expenses</center>}
      </React.Fragment>
    );
  }
}

export default Expenses;
